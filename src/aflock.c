/*
 * Implementation notes:
 *
 * - Is is assumed that all coordinates will fit into memory.
 * - Pairwise distances are computed when needed and not stored.
 * - Some parts to be parallelized.
 * - For diploid structures there are 2*nBeads beads -- ugly!
 *
 * General:
 * - Uses infinite math so DON'T compile with -ffitnite-math-only
 *
 */
#define _GNU_SOURCE

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "cf_version.h"
#include "aflock.h"
#include "ellipsoid.h"
#include "oscp.h"
#include "wio.h"


#ifdef __linux__
#include <dlfcn.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#endif


void limit_mem(size_t max_bytes)
{
#ifdef __linux__
    if(max_bytes==0)
    {
        struct sysinfo info;
        if(sysinfo(&info))
        {
            fprintf(stderr, "sysinfo returned errno: %d\n", errno);
            errno = 0;
        }
        max_bytes = info.mem_unit*info.freeram;
    }
    struct rlimit limit;

    /* The limits are not really limits... */
    getrlimit(RLIMIT_DATA, &limit);
    limit.rlim_cur = max_bytes;
    setrlimit(RLIMIT_DATA, &limit);
    return;
#else
    return;
#endif
}

#ifdef __linux__
static void* (*real_malloc)(size_t) = NULL;
static void* (*real_calloc)(size_t, size_t) = NULL;

static void mtrace_init(void)
{
    real_malloc = dlsym(RTLD_NEXT, "malloc");
    if (NULL == real_malloc) {
        fprintf(stderr, "Error in `dlsym`: %s\n", dlerror());
    }
    real_calloc = dlsym(RTLD_NEXT, "calloc");
    if (NULL == real_calloc) {
        fprintf(stderr, "Error in `dlsym`: %s\n", dlerror());
    }
    return;
}

void * malloc(size_t size)
{
    if(real_malloc==NULL) {
        mtrace_init();
    }
    void * p = real_malloc(size);
    if(p == 0)
    {
        fprintf(stderr, "aflock ran out of memory :(\n");
        fprintf(stderr, "trying to get %zu more bytes\n", size);
        exit(EXIT_FAILURE);
    }
    return p;
}

void *calloc(size_t nmemb, size_t size)
{
    if(real_calloc==NULL) {
        mtrace_init();
    }

    void *p = NULL;

    p = real_calloc(nmemb, size);
    if(p == 0)
    {
        fprintf(stderr, "aflock ran out of memory :(\n");
        fprintf(stderr, "trying to get %zu more bytes\n", nmemb*size);
        exit(EXIT_FAILURE);
    }
    return p;
}
#endif

/* Activation Distance Struct
 * To be passed to threads */
typedef struct{
    double * A; /* For storing activation distances */
    double * AD;
    size_t thread;
    size_t nThreads;
    double th_low;
    double th_high;
    fconf * fc;
    chrom * flock;
} adstruct;


/* To be passed to each thread (final_tfun)
 * when running in final mode */
typedef struct{
    size_t thread;
    size_t nThreads;
    size_t nStruct;
    size_t nBeads;
    char * wFileName;
    fconf * fc;
    chrom * flock;
    double * M; /* Counting number of captured contacts for this thread */
    double * W; /* Sum of all individual W (uint8_t) for this thread */
    double * rprof; /* Sum of radial values for this thread */
} final_tdata;


fconf * fconf_init()
{
    fconf * c = calloc(1, sizeof(fconf));
    assert(c != NULL);
    c->verbose = 1;
    c->afname = NULL;
    c->nBeads = 0;
    c->nStruct = 0;

    c->rfname = NULL;
    c->prfname = NULL;

    c->th_low = 0;
    c->th_high = 0;
    c->mode = MODE_UNKNOWN;
    c->QS = 2500;

    c->ea = -1;
    c->eb = -1;
    c->ec = -1;
    c->E = NULL;

    c->r0 = -1;
    c->vq = -1;

    c->mflock_arguments = NULL;

    c->experimental = 0;
    c->get_quality = 0;
    int nCpus = sysconf(_SC_NPROCESSORS_ONLN);
    /* Assuming this is the number of cores,
     * it is faster than using all threads */
    nCpus/= 2;
    if(nCpus < 0)
    {
        printf("WARNING: Could not figure out the number of CPUs "
               "from sysconfig!\n");
        nCpus = 4;
    }
    c->nThreads = nCpus;

    if(c->verbose > 1)
    {
        printf("Using %zu threads\n", c->nThreads);
    }

    /* Haploid or diploid cell line */
    c->diploid = 0;

    return c;
}

void fconf_free(fconf * c)
{
    if(c->E != NULL)
    {
        free(c->E);
    }
    if(c->A != NULL)
    {
        free(c->A);
    }
    if(c->afname != NULL)
    {
        free(c->afname);
    }
    if(c->rfname != NULL)
    {
        free(c->rfname);
    }
    if(c->prfname != NULL)
    {
        free(c->prfname);
    }
    if(c->mflock_arguments != NULL)
    {
        free(c->mflock_arguments);
    }
}

/* 10,000 structures, 2300 beads each, 1m 32s on MBP
 * 1m 35 s on erikwfractal. Requires loads of ram... */

static float eudist3(float * A, float * B)
{
    /* Euclidean distance between two 3D-vectors */
    return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}

static float norm3(float * X)
{
    float n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
}

void fconf_load_A(fconf * fc)
{

    fprintf(stdout, "Loading A: %s ... \n", fc->afname);

    FILE * afile = fopen(fc->afname, "r");

    if(afile == NULL)
    {
        fprintf(stderr, "Failed to open file\n");
        exit(-1);
    }

    struct stat st;
    int status = stat(fc->afname, &st);
    if(status != 0)
    {
        fprintf(stderr, "Failed to get file size of %s\n", fc->afname);
        exit(-1);
    }

    size_t fsize = st.st_size;

    fc->A = malloc(fsize);
    if(fc->A == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for A\n");
        exit(EXIT_FAILURE);
    }

    size_t nRead = fread(fc->A, sizeof(double), fsize/sizeof(double), afile);

    if(nRead*sizeof(double) != fsize)
    {
        fprintf(stderr, "nRead = %zu, fsize = %zu\n", nRead, fsize);
        exit(-1);
    }

    fclose(afile);

    for(size_t kk = 0 ; kk<fsize/sizeof(double) ; kk++)
    {
        if(isfinite(fc->A[kk]))
        {
            if(fc->A[kk] < 0)
            {
                fprintf(stderr, "ERROR: A[%zu] = %f < 0\n", kk, fc->A[kk]);
                exit(-1);
            }
            if(fc->A[kk] > 1)
            {
                fprintf(stderr, "ERROR: A[%zu] = %f > 1\n", kk, fc->A[kk]);
                exit(-1);
            }
        }
    }


    fc->nBeads = sqrt(fsize/sizeof(double));
    assert(pow(fc->nBeads, 2) == fsize/sizeof(double));
    fprintf(stdout, "   .. contains [%zu x %zu] elements.\n", fc->nBeads, fc->nBeads);
    return;
}

int chrom_init(chrom * c, size_t nQ, size_t n)
{
    c->X = NULL;
    c->nQ = 0; /* Number of elements in the queue */
    c->Q = malloc(2*nQ*sizeof(uint32_t));
    assert(c->Q != NULL);

    size_t bsize = 64;
    c->wfName = malloc(bsize);
    assert(c->wfName != NULL);
    snprintf(c->wfName, bsize, "cf_%06zu/W.uint8.gz", n+1);
    c->W = 0;

    if(c->wfName == NULL)
        return 1;

    c->xfName = malloc(bsize);
    assert(c->xfName != NULL);
    snprintf(c->xfName, bsize, "cf_%06zu/coords.csv", n+1);

    if(c->xfName == NULL)
        return 1;

    return 0;
}

void ch_free(chrom * c)
{
    if(c->W != NULL)
        free(c->W);
    if(c->X != NULL)
        free(c->X);
    if(c->Q != NULL)
        free(c->Q);
    free(c->wfName);
    free(c->xfName);
    c->nQ = 0;
    return;
}


int ch_load_X(fconf * fc, chrom * c)
{

    size_t nBeads = fc->nBeads*(1+fc->diploid);

    c->X = malloc(nBeads*3*sizeof(float));
    assert(c->X != NULL);
    /* Read from c->xfName ... */

    //  fprintf(stdout, "Reading X-data from %s\n", c->xfName);

    FILE * f = fopen(c->xfName, "r");
    if(f == NULL)
    {
        fprintf(stderr, "\rCan't open %s\n", c->xfName);
        fprintf(stderr, "Did you forget to run mflock after aflock -I?\n");
        exit(-1);
    }

    char * line = malloc(1024*sizeof(char));
    assert(line != NULL);
    size_t len = 1024*sizeof(char);

    char delim[] = ",";
    for(size_t ll = 0; ll < nBeads; ll++)
    {
        int read = getline(&line, &len, f);
        if(read == -1)
        {
            printf("Failed to read line %zu\n", ll+1);
            exit(-1);
        }
        char *ptr = strtok(line, delim);
        c->X[3*ll] = atof(ptr);
        ptr = strtok(NULL, delim);
        c->X[3*ll+1] = atof(ptr);
        ptr = strtok(NULL, delim);
        c->X[3*ll+2] = atof(ptr);
    }

    free(line);
    fclose(f);

    return 0;
}

int cmp_float_reverse(const void * A, const void * B)
{
    float * dA = (float *) A;
    float * dB = (float *) B;

    if(dA[0] > dB[0])
        return 1;
    if(dA[0] < dB[0])
        return -1;
    return 0;
}

int cmp_float(const void * A, const void * B)
{
    float * dA = (float *) A;
    float * dB = (float *) B;

    if(dA[0] > dB[0])
        return 1;
    if(dA[0] < dB[0])
        return -1;
    return 0;
}

void chrom_assign(fconf * fc, chrom * ch, size_t kk, size_t ll)
{
    if(ch->nQ == fc->QS)
    {
        /* Empty the queue first by writing to disk. */
        struct_write_W(fc, ch);
    }

    /* Insert in Queue */
    ch->Q[2*ch->nQ] = kk;
    ch->Q[2*ch->nQ+1] = ll;

    ch->nQ++;
}


chrom * load_structures(fconf * fc)
{
    fprintf(stdout, "Initializing %zu structures ...\n", fc->nStruct);
    chrom * flock = malloc(fc->nStruct*sizeof(chrom));
    assert(flock != NULL);

    /* Load them */
    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
        printf("\r%zu", kk+1); fflush(stdout);
        chrom_init(&flock[kk], fc->QS, kk);
    }
    printf("\r");

    if(fc->mode == MODE_INIT)
    {
        printf("No coordinates to be loaded\n");
    }
    else
    {
        printf("Loading coordinates ...\n");

        for(size_t kk = 0; kk<fc->nStruct; kk++)
        {
            printf("\r%zu", kk); fflush(stdout);
            ch_load_X(fc, &flock[kk]);
        }
        printf("\r");
    }

    printf("\r");

    fprintf(stdout, "   ... done.\n");
    return flock;
}

void writew(char * wfname, uint8_t * W, size_t nElements)
{
    wio_write(wfname, nElements, W);
}

void struct_write_W0(fconf * fc, chrom * cc, uint8_t * W0)
{
    if(fc->diploid == 1)
    {
        wio_write(cc->wfName, pow(2*fc->nBeads,2), W0);
    } else {
        wio_write(cc->wfName, pow(fc->nBeads,2), W0);
    }
}

void struct_write_R0(fconf * fc, chrom * cc, double * R0)
{
    if(fc->diploid == 1)
    {
        wio_write(cc->rfName, 2*fc->nBeads*sizeof(double), (void *) R0);
    } else {
        wio_write(cc->rfName, fc->nBeads*sizeof(double), (void *) R0);
    }
}

uint8_t * readw(char * wfName, size_t nElements)
{
    return wio_read(wfName, &nElements);
}

void struct_write_W_AD(fconf * fc,
                       chrom * c,
                       double * AD,
                       double th_high,
                       double th_low)
{
    /* Write a new W matrix to disk, setting contacts between beads to
     * 1 if the distance is below the activation distance given by AD
     *
     * This will read an existing W matrix and only update contacts for which
     * A is in the threshold range.
     * */


    size_t nBeads = (1+fc->diploid)*fc->nBeads;
    size_t N = fc->nBeads;

    size_t nRead = pow(nBeads, 2);
    uint8_t * W = readw(c->wfName, nRead);
    double * A = fc->A;

    for(size_t aa = 0; aa < nBeads; aa++)
    {
        for(size_t bb = aa+1; bb < nBeads; bb++)
        {
            size_t idx1 = (aa % N) + (bb % N)*fc->nBeads;
            //size_t idx2 = (bb % N) + (aa % N)*fc->nBeads;
            size_t widx1 = aa + bb*nBeads;
            size_t widx2 = bb + aa*nBeads;

            /* If within the specified theta range */
            if( (A[idx1] < th_high) & (A[idx1] >= th_low) )
            {
                /* And a activation distance is specified */
                if(isfinite(AD[idx1]))
                {
                    /* Write 1 to W if the bead pair distance
                       for this structure is within the activation distance */
                    double d = eudist3(c->X+aa*3, c->X+bb*3);

                    if(d<AD[idx1])
                    {
                        W[widx1] = 1;
                        W[widx2] = 1;
                    } else {
                        W[widx1] = 0;
                        W[widx2] = 0;
                    }
                } else {
                    W[widx1] = 0;
                    W[widx2] = 0;
                }
            }
        }
    }
    /* Write back to disk */
    //  printf("Will write to: %s\n", c->wfName); fflush(stdout);
    wio_write(c->wfName, nRead, W);
    free(W);
    return;
}

void struct_write_W(fconf * fc, chrom * c)
{
    //  fprintf(stdout, "Will open %s and write %zu new constraints\n", c->wfName, c->nQ);

    size_t nRead = pow(fc->nBeads, 2);
    uint8_t * W = readw(c->wfName, nRead);


    /* Update W with new contacts */
    for(size_t pp = 0; pp<c->nQ; pp++)
    {
        size_t kk = c->Q[2*pp];
        size_t ll = c->Q[2*pp+1];

        W[kk+fc->nBeads*ll] = 1;
        W[ll+fc->nBeads*kk] = 1;
    }

    /* Write back to disk */
    wio_write(c->wfName, nRead, W);
    c->nQ = 0;

    free(W);

    return;
}

void flock_reset(fconf * fc, chrom * flock, double th_high, double th_low)
{
    /* Reset all contacts in each W where th_high < A <= th_low
     */
    fprintf(stdout, "Resetting contacts where %f <= A < %f\n", th_low, th_high);
    size_t nElements = (size_t) pow(fc->nBeads, 2);

    for(size_t pp = 0; pp<fc->nStruct; pp++)
    {
        printf("\r%zu", pp); fflush(stdout);

        uint8_t * W = readw(flock[pp].wfName, nElements);

        for(size_t kk = 0; kk<nElements; kk++)
        {
            if(fc->A[kk] >= th_low && fc->A[kk] < th_high)
            {
                W[kk] = 0;
            }
        }

        writew(flock[pp].wfName, W, nElements);

        free(W);
    }
    printf("\r");
    printf("Done!  \n");


}

void flock_assign_with_w_randblock(fconf * fc, chrom * flock, size_t kk, size_t ll, float prob)
{
    /* Assign a contacts between bead kk and ll with with probability prob
     */

    size_t nChrom = fc->nStruct;
    size_t nAssign = round(prob*(float) nChrom);

    /* Array of pairwise distance between point kk and ll */
    float * D = malloc(nChrom*sizeof(float));
    assert(D != NULL);

    for(size_t pp = 0; pp<nChrom; pp++)
    {
        D[pp] = eudist3(flock[pp].X+kk*3, flock[pp].X+ll*3);
        if(flock[pp].W[kk+fc->nBeads*ll] == 1)
        {
            double drand = (double) rand()/ (double) RAND_MAX;
            if(D[pp] > 2.5*fc->r0)
            {
                D[pp] = drand + 2;
            }
            else
            {
                D[pp] = drand;
            }

        }
    }

    float * DS = malloc(nChrom*sizeof(float));
    assert(DS != NULL);
    memcpy(DS, D, nChrom*sizeof(float));

    qsort(DS, nChrom, sizeof(float), cmp_float);

    float threshold = DS[nAssign-1];
    free(DS);

    if(0)
    {
        fprintf(stdout, "Found threshold: %f\n", threshold);
    }

    size_t nAssigned = 0;

    for(size_t pp = 0; pp<nChrom; pp++)
    {
        if((D[pp]<=threshold) && (nAssigned < nAssign))
            /* To avoid the case where multiple structures with the same distance */
        {
            flock[pp].W[kk+fc->nBeads*ll] = 1;
            flock[pp].W[ll+fc->nBeads*kk] = 1;
            nAssigned++;
        } else {
            flock[pp].W[kk+fc->nBeads*ll] = 0;
            flock[pp].W[ll+fc->nBeads*kk] = 0;
        }
    }

    if(nAssigned != nAssign)
    {
        fprintf(stderr, "To assign: %zu\n", nAssign);
        fprintf(stderr, "Assigned: %zu\n", nAssigned);
        fprintf(stderr, "Threshold: %f\n", threshold);
        assert(0);
    }
    free(D);
    return;
}

void flock_assign_with_w(fconf * fc, chrom * flock, size_t kk, size_t ll, float prob)
{
    /* Assign a contacts between bead kk and ll with with probability prob
     */

    size_t nChrom = fc->nStruct;
    size_t nAssign = round(prob*(float) nChrom);

    /* Array of pairwise distance between point kk and ll */
    float * D = malloc(nChrom*sizeof(float));
    assert(D != NULL);
    for(size_t pp = 0; pp<nChrom; pp++)
    {
        D[pp] = eudist3(flock[pp].X+kk*3, flock[pp].X+ll*3);
        if(flock[pp].W[kk+fc->nBeads*ll] == 1)
        {
            if(D[pp] > 2.5*fc->r0)
            {
                D[pp] += 2;
            }
        }
    }

    float * DS = malloc(nChrom*sizeof(float));
    assert(DS != NULL);
    memcpy(DS, D, nChrom*sizeof(float));

    qsort(DS, nChrom, sizeof(float), cmp_float);

    float threshold = DS[nAssign-1];
    free(DS);

    if(0)
    {
        fprintf(stdout, "Found threshold: %f\n", threshold);
    }

    size_t nAssigned = 0;

    for(size_t pp = 0; pp<nChrom; pp++)
    {
        if((D[pp]<=threshold) && (nAssigned < nAssign))
            /* To avoid the case where multiple structures with the same distance */
        {
            flock[pp].W[kk+fc->nBeads*ll] = 1;
            flock[pp].W[ll+fc->nBeads*kk] = 1;
            nAssigned++;
        } else {
            flock[pp].W[kk+fc->nBeads*ll] = 0;
            flock[pp].W[ll+fc->nBeads*kk] = 0;
        }
    }

    if(nAssigned != nAssign)
    {
        fprintf(stderr, "To assign: %zu\n", nAssign);
        fprintf(stderr, "Assigned: %zu\n", nAssigned);
        fprintf(stderr, "Threshold: %f\n", threshold);
        assert(0);
    }
    free(D);
    return;
}

void flock_assign(fconf * fc, chrom * flock, size_t kk, size_t ll, float prob)
{
    /* Assign a contacts between bead kk and ll with with probability prob
     */

    size_t nChrom = fc->nStruct;
    size_t nAssign = round(prob*(float) nChrom);
    assert(nAssign > 0);
    assert(nAssign <= nChrom);
    if(nAssign == 0)
    {
        return;
    }

    if(0)
    {
        printf("prob=%f -> %zu/%zu structures\n", prob, nAssign, nChrom);
    }
    float * D = malloc(nChrom*sizeof(float));
    assert(D != NULL);
    for(size_t pp = 0; pp<nChrom; pp++)
    {
        D[pp] = eudist3(flock[pp].X+kk*3, flock[pp].X+ll*3);
        /* If we knew which contacts that failed, we could increase the
         * distance of them here in order not to try them again...
         * (at least not in the next round). */
    }

    float * DS = malloc(nChrom*sizeof(float));
    assert(DS != NULL);
    memcpy(DS, D, nChrom*sizeof(float));

    qsort(DS, nChrom, sizeof(float), cmp_float);


    float threshold = DS[nAssign-1];

    free(DS);

    if(0)
    {
        fprintf(stdout, "Found threshold: %f\n", threshold);
    }

    size_t nAssigned = 0;

    /* (1) In case that the distances are equal for several structures
     * the contact is assigned to the first of the structures
     */

    if(fc->experimental == 0)
    {
        for(size_t pp = 0; pp<nChrom; pp++)
        {
            if((D[pp] <= threshold) &&
               (nAssigned < nAssign)) // (1)
            {
                chrom_assign(fc, flock+pp, kk, ll);
                nAssigned++;
            }
        }
    } else {
        /* Experimental random assignment! */

        uint8_t * set = malloc(nChrom*sizeof(uint8_t));
        assert(set != NULL);
        memset(set, 0, nChrom*sizeof(uint8_t));

        for(size_t pp =0; pp<nAssign; pp++)
        {
            int isSet = 0;
            while(isSet == 0)
            {
                size_t pos = rand() % nChrom;
                if(set[pos] == 0)
                {
                    isSet = 1;
                    set[pos] = 1;
                    chrom_assign(fc, flock+pos, kk, ll);
                    nAssigned++;
                }
            }
        }
        free(set);
    }
    assert(nAssigned == nAssign);
    free(D);

    return;
}


static void usage()
{
    printf("Accepted arguments:\n");
    printf(" -A file.dat, --Afile file.dat\n\tcontact probability matrix, square double.\n");
    printf(" -n nStruct, --nStruct nStruct\n\tspecify the number of structures\n");
    printf(" -h th, --high th\n\thigher threshold\n");
    printf(" -l th, --low th\n\tlower threshold\n");
    printf(" -R radius, --radius radius\n\tgive bead radius (-Q preferred)\n");
    printf(" -Q vq, --vq vq\n\tgive volume quotient of beads/nuclei (instead of -R)\n");
    printf(" --ea a --eb b --ec c\n\tSpecify axes of ellipsoidal (otherwhise sphere). 1=a>=b>=c>0\n");
    printf(" -P args, --mArgs args\n\tcommand line arguments to pass to mflock\n");
    printf(" --rpos RPOS.double\n\tsupply radial preference for each bead\n");
    printf(" --prpos PRPOS.double\n\tsupply the probability that each bead is given the radial constraint\n");
    printf(" --diploid\n\tDiploid cell line (otherwise haploid)\n");
    printf(" --threads #\n\tSet number of threads to use\n");
    printf("MODES:\n");
    printf(" -I, --init\n\tinitialize new structures\n");
    printf(" -F, --final\n\tfinal stuff\n");
    printf(" -U, --update\n\tupdate\n");
    printf(" -B, --blockupdate\n\tupdate and block failed contacts\n");
    printf(" -X, --update_experimental\n\texperimental update algorithm\n");
    printf(" -S, --score\n\tscore current quality (not for -I)\n");
    printf("See the man page for more information\n");

}

static int argparsing(fconf * p, int argc, char ** argv)
{

    struct option longopts[] = {
        // MODES
        { "help",         no_argument,       NULL,  'u' },
        { "init",         no_argument,       NULL,  'I' },
        { "update",       no_argument,       NULL,  'U' },
        { "blockupdate",  no_argument,       NULL,  'B' },
        { "final",        no_argument,       NULL,  'F' },
        { "score",        no_argument,       NULL,  'S' },
        { "update_experimental",
          no_argument,       NULL,  'X' },
        /* Data */
        { "Afile",        required_argument, NULL,  'A' },
        { "Rfile",        required_argument, NULL,  'R' },
        /* Settings */
        { "mArgs",        required_argument, NULL,  'P' },
        { "nStruct",      required_argument, NULL,  'n' },
        { "high",         required_argument, NULL,  'h' },
        { "low",          required_argument, NULL,  'l' },
        { "threads",      required_argument, NULL,  'T' },
        /* Geometry */
        { "ea",           required_argument, NULL,  'a'},
        { "eb",           required_argument, NULL,  'b'},
        { "ec",           required_argument, NULL,  'c'},
        { "radius",       required_argument, NULL,  'R' },
        { "vq",           required_argument, NULL,  'Q' },
        {"rpos",          required_argument, NULL,  'r' },
        {"prpos",         required_argument, NULL,  'p' },
        {"version",       no_argument,       NULL,  'i' },
        {"diploid",       no_argument,       NULL,  'D' },
        { NULL,           0,                 NULL,   0  }
    };

    int ch;
    while((ch = getopt_long(argc, argv,
                            "A:BDEFIP:Q:R:SUXh:l:n:p:r:u",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'i':
            printf("aflock (chromflock version %s)\n", cf_version);
            printf("Build date: %s, %s\n", __DATE__, __TIME__);
            printf("GIT HASH: %s\n", GIT_VERSION);
            printf("Compiler: %s\n", CC_VERSION);
            exit(0);
        case 'A': /* sphere contact probability matrix */
            if(p->afname != NULL) { free(p->afname); }
            p->afname = malloc(strlen(optarg)+1);
            assert(p->afname != NULL);
            strcpy(p->afname, optarg);
            break;
        case 'B':
            p->mode = MODE_UPDATE_BLOCKED;
            break;
        case 'D':
            p->diploid = 1;
            break;
        case 'E':
            p->experimental = 1;
            break;
        case 'F':
            p->mode = MODE_FINAL;
            break;
        case 'I': /* init, i.e., create W */
            p->mode = MODE_INIT;
            break;
        case 'P':
            if(p->mflock_arguments != NULL) { free(p->mflock_arguments); }
            p->mflock_arguments = malloc(strlen(optarg)+1);
            assert(p->mflock_arguments != NULL);
            strcpy(p->mflock_arguments, optarg);
            break;
        case 'Q':
            p->vq = atof(optarg);
            break;
        case 'R':
            p->r0 = atof(optarg);
            break;
        case 'S':
            p->get_quality = 1;
            break;
        case 'T':
            p->nThreads = atoi(optarg);
            printf("Using %zu threads\n", p->nThreads);
            break;

        case 'U':
            p->mode = MODE_UPDATE;
            break;
        case 'X':
            p->mode = MODE_EXPERIMENTAL;
            break;
        case 'a':
            p->ea = atof(optarg);
            break;
        case 'b':
            p->eb = atof(optarg);
            break;
        case 'c':
            p->ec = atof(optarg);
            break;
        case 'h': /* high threshold */
            p->th_high = atof(optarg);
            break;
        case 'l': /* low threshold */
            p->th_low = atof(optarg);
            break;
        case 'n':
            p->nStruct = atol(optarg);
            break;
        case 'p':
            if(p->prfname != NULL) { free(p->prfname); }
            p->prfname = malloc(strlen(optarg)+1);
            assert(p->prfname != NULL);
            strcpy(p->prfname, optarg);
            break;
        case 'r':
            if(p->rfname != NULL) { free(p->rfname); }
            p->rfname = malloc(strlen(optarg)+1);
            assert(p->rfname != NULL);
            strcpy(p->rfname, optarg);
            break;
        case 'u':
            usage();
            exit(0);
        default:
            return(1);
        }
    }

    if(p->afname == NULL)
    {
        printf("No A file was specified\n");
        return 1;
    }

    if(p->nStruct == 0)
    {
        printf("Don't know how many structures to deal with\n");
        return 1;
    }

    if((p->rfname != NULL) & (p->prfname == NULL))
    {
        printf("--rpos given but not --prpos\n");
        return 1;
    }

    if((p->rfname == NULL) & (p->prfname != NULL))
    {
        printf("--prpos given but not --rpos\n");
        return 1;
    }

    return 0;
}


void struct_get_interactions(fconf * fc, chrom * cc, double * M,
                             double cDistance)
/* cDistance: Contact distance in normalized units */
{
    size_t N = (1+fc->diploid)*fc->nBeads;

    for(size_t kk = 0; kk<N; kk++)
    {
        for(size_t ll = kk+1; ll<N; ll++)
        {
            if( eudist3(cc->X + 3*kk, cc->X + 3*ll) <= cDistance )
            {
                M[kk + ll*N]++;
                M[ll + kk*N]++;
            }
        }
    }
    return;
}


void struct_add_radius(fconf * fc, chrom * cc, double * R)
{
    /* Calculate radius for each point in structure */

    size_t nBeads = (1+fc->diploid)*fc->nBeads;

    /* For spherical geometry */
    if(fc->E == NULL)   {
        for(size_t pp = 0; pp<nBeads; pp++)
        {
            R[pp]+=norm3(cc->X+3*pp);
        }
    }

    /* For ellipsoidal geometry */
    if(fc->E != NULL) {
        double XT[3];
        for(size_t pp = 0; pp<nBeads; pp++)
        {
            XT[0] = (double) cc->X[3*pp + 0];
            XT[1] = (double) cc->X[3*pp + 1];
            XT[2] = (double) cc->X[3*pp + 2];

            R[pp]+=elli_getScale(fc->E, XT);
        }
    }
    return;
}


void set_bead_radius(fconf * fc)
{
    double captureFactor = 4; /* TODO: deduce from dynamics program  */

    if((fc->vq > 0) && (fc->r0 > 0)) {
        fprintf(stderr, "Either -R or -Q has to be specified, not both!\n");
        exit(-1);
    }

    size_t nBeads = fc->nBeads;
    if(fc->diploid == 1)
    {   nBeads *= 2;    }

    double vDomain = 4.0/3.0*M_PI;
    if(fc->E != NULL)
    {
        vDomain = elli_vol(fc->E);
    }

    int allSet = 0;

    if(fc->vq > 0) {
        printf("Setting beads radius from volume quotient\n");
        if( (fc->vq <= 0) || (fc->vq > 1)) {
            fprintf(stderr, "The volume quotient has be be in ]0, 1], use .4 if unsure\n");
            exit(-1);
        }
        fc->r0 = cbrt(fc->vq*vDomain / (nBeads*4.0*M_PI/3.0));
        fc->dContact = captureFactor*fc->r0;
        allSet = 1;
    }

    if(fc->r0 > 0 && allSet == 0) {
        printf("Setting volume quotient from bead radius\n");
        if(fc->r0 > 1) {
            fprintf(stderr, "A radius of %f does not make sense\n", fc->r0);
            exit(-1);
        }
        allSet = 1;
    }

    assert(allSet == 1);
    double vBeads = (double) nBeads*pow(fc->r0, 3)*4.0/3.0*M_PI;
    fc->vq = vBeads/vDomain;
    fc->dContact = captureFactor*fc->r0;
    fprintf(stdout, "   Number of beads: %d*%zu = %zu\n", fc->diploid + 1, fc->nBeads, nBeads);
    fprintf(stdout, "   Bead radius: %f\n", fc->r0);
    fprintf(stdout, "   Total bead volume: %f\n", vBeads);
    fprintf(stdout, "   Sphere volume: %f\n", vDomain);
    fprintf(stdout, "   Volume Quotient: %f\n", vBeads/vDomain);
    fprintf(stdout, "   Contact distance (%.1f x bead radius): %f\n", captureFactor, fc->dContact);
    assert(fabs(fc->vq - vBeads/vDomain)<1e-6);
    assert(vBeads/vDomain < 1);
    return;
}

double structQuality(fconf * fc, chrom * cc)
{

    size_t W0S0 = 0;
    size_t W0S1 = 0;
    size_t W1S0 = 0;
    size_t W1S1 = 0;

    /* load W */
    uint8_t * W = readw(cc->wfName, pow(fc->nBeads, 2));

    size_t N = fc->nBeads;

    for(size_t kk = 0; kk<N; kk++)
    {
        for(size_t ll = kk+1; ll<N; ll++)
        {
            int contact = 0; /* if the beads are in contact */
            int set = 0; /* if there should be a contact */

            double d = eudist3(cc->X + 3*kk, cc->X + 3*ll);

            if(d < fc->dContact)
            {
                contact = 1;
            }

            if(W[kk*N+ll] == 1)
            {
                set = 1;
            }

            if(set == 1)
            {
                if(contact == 1)
                {
                    W1S1++;
                } else {
                    W1S0++;
                }
            } else {
                if(contact == 1)
                {
                    W0S1++;
                } else
                {
                    W0S0++;
                }
            }
        }

    }

    free(W);

    printf("W==0,S==0, %zu, W==0,S==1, %zu, W==1, S==0, %zu, W==1, S==1 %zu\n", W0S0, W0S1, W1S0, W1S1);

    //return ((double) W0S0 + (double) W1S1)/((double) W0S0 + (double) W0S1 + (double) W1S0 + (double) W1S1);
    return (double) W1S1 / ( (double) W1S0 + (double) W1S1 ); // got / wanted
}

double getQuality(fconf * fc, chrom * flock)
{

    double * quality = malloc(fc->nStruct*sizeof(double));
    assert(quality != NULL);
    memset(quality, 0, fc->nStruct*sizeof(double));

    double q = 0;

    for(size_t cc = 0; cc < fc->nStruct; cc++)
    {
        quality[cc] = structQuality(fc, flock+cc);
        q+=quality[cc];
    }

    free(quality);
    return q / fc->nStruct;
}

void echo_args(int argc, char ** argv)
{
    for(int kk = 0; kk<argc; kk++)
    {
        printf("%s ", argv[kk]);
    }
    printf("\n");
    return;
}

void calc_activation_distance_th( adstruct * s)
{
    fconf * fc = s->fc;
    double * A = s->A;
    double * AD = s->AD;
    size_t thread = s->thread;
    size_t nThreads = s->nThreads;
    double th_low = s->th_low;
    double th_high = s->th_high;
    chrom * flock = s->flock;

    size_t nDist = (1 + 3*fc->diploid)*fc->nStruct;
    size_t B = fc->nBeads;

    /* temporary array for each bead pair */
    float * DS = malloc(nDist*sizeof(float));
    assert(DS != NULL);

    for(size_t aa = thread; aa < fc->nBeads; aa = aa+nThreads)
    {

        for(size_t bb = aa+1; bb < fc->nBeads; bb++)
        {
            size_t idx1 = aa+bb*fc->nBeads;
            size_t idx2 = bb+aa*fc->nBeads;

            if((A[idx1] < th_high) & (A[idx1] >= th_low))
            {
                /* Only difference between diploid and haploid is that we want twice
                 * as many contacts for the diploid */
                size_t nAssign = (1+fc->diploid)*round(A[idx1]*(float) fc->nStruct);

                if(nAssign>0)
                {
                    /* Get distance between bead aa and bb in all structures */
                    if(fc->diploid == 0) {
                        for(size_t ss = 0; ss< fc->nStruct; ss++) {
                            DS[ss] = eudist3(flock[ss].X+aa*3, flock[ss].X+bb*3);
                        }}

                    if(fc->diploid == 1) {
                        for(size_t ss = 0; ss< fc->nStruct; ss++) {
                            DS[4*ss+0] = eudist3(flock[ss].X + (0+aa)*3, flock[ss].X + (0+bb)*3);
                            DS[4*ss+1] = eudist3(flock[ss].X + (B+aa)*3, flock[ss].X + (0+bb)*3);
                            DS[4*ss+2] = eudist3(flock[ss].X + (0+aa)*3, flock[ss].X + (B+bb)*3);
                            DS[4*ss+3] = eudist3(flock[ss].X + (B+aa)*3, flock[ss].X + (B+bb)*3);
                        }}

                    qsort(DS, nDist, sizeof(float), cmp_float);
                    double act_dist = 10;
                    if(nAssign < nDist)
                    {
                        act_dist = (DS[nAssign-1] + DS[nAssign])/2;
                    }
                    AD[idx1] = act_dist;
                    AD[idx2] = act_dist;
                }
            }
        }
    }

    free(DS);
}


void * final_tfun(void * data)
{
    final_tdata * td = (final_tdata *) data;
    const size_t N = td->nBeads;

    /* This thread sums up all individual W to td->W */
    td->W = malloc(pow(N, 2)*sizeof(double));
    assert(td->W != NULL);

    if(td->thread == 0){
        printf("   Summing up wanted contacts ... \n");
    }
    td->wFileName = malloc(1024*sizeof(char));
    assert(td->wFileName != NULL);

    for(size_t pp = td->thread; pp < td->nStruct ; pp+=td->nThreads)
    {
        sprintf(td->wFileName, "cf_%06zu/W.uint8.gz", pp+1); fflush(stdout);
        size_t nRead = 0;
        uint8_t * W = (uint8_t *) wio_read(td->wFileName, &nRead);
        if(nRead != pow(N, 2))
        {
            printf("Error: Read %zu bytes from %s, expected %zu\n",
                   nRead,
                   td->wFileName,
                   (size_t) pow(td->nBeads, 2));
            exit(1);
        }
        for(size_t kk = 0; kk<N*N; kk++)
        {
            td->W[kk] += W[kk];
        }
        free(W);
    }

    if(td->thread == 0)
    {
        printf("   Finding beads in contact ... \n");
    }
    /* TODO: dContact should be deduced from the dynamics program */
    const double cDistance = td->fc->dContact;
    for(size_t pp = td->thread; pp < td->nStruct ; pp+=td->nThreads)
    {
        float * X = td->flock[pp].X;
        for(size_t kk = 0; kk<N; kk++)
        {
            for(size_t ll = kk+1; ll<N; ll++)
            {
                if( eudist3(X + 3*kk, X + 3*ll) <= cDistance )
                {
                    td->M[kk + ll*N]++;
                }
            }
        }
    }

    if(td->thread == 0)
    {
        printf("   Constructing radial profile ... \n");
    }
    td->rprof = calloc(N, sizeof(double));
    assert(td->rprof != NULL);
    for(size_t pp = td->thread; pp < td->nStruct ; pp+=td->nThreads)
    {
        float * X = td->flock[pp].X;
        if(td->fc->E == NULL)
        {
            for(size_t kk = 0; kk<N; kk++)
            {
                td->rprof[kk]+=norm3(X+3*kk);
            }
        }

        /* For ellipsoidal geometry */
        if(td->fc->E != NULL) {
            double XT[3];
            for(size_t pp = 0; pp<N; pp++)
            {
                XT[0] = (double) X[3*pp + 0];
                XT[1] = (double) X[3*pp + 1];
                XT[2] = (double) X[3*pp + 2];
                td->rprof[pp]+=elli_getScale(td->fc->E, XT);
            }
        }
    }

    free(td->wFileName);

    return NULL;
}


void * run_struct_write_W_AD(void * P)
{
    adstruct * s = (adstruct *) P;
    for(size_t kk = s->thread; kk < s->fc->nStruct ; kk += s->nThreads)
    {
        struct_write_W_AD(s->fc, &s->flock[kk], s->AD, s->th_high, s->th_low);
    }

    return NULL;
}


void * run_calc_activation_distance_th(void * P)
{
    fflush(stdout);
    adstruct * s = (adstruct *) P;
    calc_activation_distance_th(s);
    return NULL;
}


void calc_activation_distance(fconf * fc, chrom * flock, double th_high, double th_low, double * AD)
{
    /* Consider all positions in PR where p<th_high and p>= th_low
     * For each beads:
     * hand the constrains to the structures that has the closest radii already
     * similar to how the contact restraints are handled
     */

    printf("%s\n", __func__);
    double * A = fc->A;

    /* Each thread will process elements thread+kk*nThreads, kk = 0, 1, ...
     */

    pthread_t * threads = malloc(fc->nThreads*sizeof(pthread_t));
    assert(threads != NULL);
    adstruct ** S = malloc(fc->nThreads*sizeof(adstruct * ));
    assert(S != NULL);
    for(size_t kk = 0; kk<fc->nThreads; kk++)
    {
        S[kk] = malloc(sizeof(adstruct));
        assert(S[kk] != NULL);
        S[kk]->thread = kk;
        S[kk]->nThreads = fc->nThreads;

        S[kk]->AD = AD;
        S[kk]->A = A;
        S[kk] -> th_high = th_high;
        S[kk] -> th_low = th_low;
        S[kk]->flock = flock;
        S[kk]->fc = fc;

        pthread_create(&threads[kk],
                       NULL,
                       run_calc_activation_distance_th,
                       (void *) S[kk]);
    }

    for(size_t kk = 0; kk<fc->nThreads; kk++)
    {
        pthread_join(threads[kk], NULL);
        free(S[kk]);
    }

    free(S);
    free(threads);

    fprintf(stdout, "\n\n");
    fprintf(stdout, "Writing activation distances to 'activationDistance.double'\n");
    wio_write("activationDistance.double",
              fc->nBeads*fc->nBeads*sizeof(double),
              (void *) AD);
    fprintf(stdout, "Done\n"); fflush(stdout);

    return;
}

void flock_updateR(fconf * fc, chrom * flock, double th_high, double th_low)
{
    /* Consider all positions in PR where p<th_high and p>= th_low
     * For each beads:
     * hand the constrains to the structures that has the closest radii already
     * similar to how the contact restraints are handled
     */

    /* To start with, all coordinates should already be available in flock. */

    /* Load preferred radial positions R
     * Load probability/proportion of structures having the constraint into PR */

    printf("%s\n", __func__);

    /* Read R: radius per bead and PR: probability of use of each bead */
    size_t nbytes_r = 0;
    size_t nbytes_pr = 0;
    double * R = (double *) wio_read(fc->rfname, &nbytes_r);
    double * PR = (double *) wio_read(fc->prfname, &nbytes_pr);

    assert(nbytes_r == nbytes_pr);
    assert(nbytes_r == fc->nBeads*sizeof(double));

    /* Find a threshold for each positions */
    double * TH = (double * ) malloc(fc->nBeads*sizeof(double));
    assert(TH != NULL);
    float * DS = malloc(fc->nStruct*sizeof(float));
    assert(DS != NULL);
    for(size_t pp = 0; pp < fc->nBeads; pp++)
    {
        if((PR[pp] < th_high) & (PR[pp] >= th_low))
        {
            size_t nAssign = round(PR[pp]*(float) fc->nStruct);
            if(nAssign>0)
            {

                for(size_t ss = 0; ss< fc->nStruct; ss++)
                {
                    DS[ss] = norm3(flock[ss].X+pp*3);
                    //printf("%f ", DS[ss]);
                }
                //printf("\n");
                qsort(DS, fc->nStruct, sizeof(float), cmp_float);
                TH[pp] = DS[nAssign-1];
                if(0){
                    printf("PR[%zu]=%f, R[%zu]=%f. nassign: %zu. Th[%zu]=%f\n",
                           pp, PR[pp], pp, R[pp], nAssign, pp, TH[pp]);
                }
            }
        }
    }
    /* Read/Write all individual R-files. */
    for(size_t ss = 0; ss < fc->nStruct ; ss++)
    {
        flock[ss].rfName = malloc(1024*sizeof(char));
        assert(flock[ss].rfName != NULL);
        sprintf(flock[ss].rfName, "cf_%06zu/radius.double.gz", ss+1);

        size_t bytesRead = 0;
        double * SR = (double *) wio_read( flock[ss].rfName, &bytesRead );
        if(bytesRead != sizeof(double)*fc->nBeads)
        {
            fprintf(stderr, "%d ERROR: Read %zu bytes from %s, expected %zu.\n", __LINE__,
                    bytesRead, flock[ss].rfName, fc->nBeads*sizeof(double));
            exit(-1);
        }

        /* Go through all pp "beads" of structure ss */
        for(size_t pp = 0 ; pp < fc->nBeads ; pp++)
        {
            if((PR[pp] < th_high) & (PR[pp] >= th_low))
                /* Don't change the values for the other beads */
            {
                double r = norm3(flock[ss].X+pp*3);
                if(r <= TH[pp])
                {
                    SR[pp] = R[pp];
                } else {
                    SR[pp] = NAN;
                }
            }
        }
        wio_write(flock[ss].rfName,
                  fc->nBeads*sizeof(double),
                  (void *) SR);
        free(SR);
        free(flock[ss].rfName);
    }

    free(DS);
    free(TH);
    free(PR);
    free(R);
    return;
}

uint8_t * initial_W(fconf * fc)
{
    if(fc->diploid == 0)
    {
        uint8_t * W0 = malloc(pow(fc->nBeads, 2)*sizeof(uint8_t));
        assert(W0 != NULL);
        for(size_t kk = 0; kk<pow(fc->nBeads,2); kk++)
        {
            if(fc->A[kk] == 1)
            {
                W0[kk] = 1;
            } else {
                W0[kk] = 0;
            }
        }
        return W0;
    }

    if(fc->diploid == 1)
    {
        uint8_t * W0 = calloc(pow(2*fc->nBeads, 2), sizeof(uint8_t));
        assert(W0 != NULL);
        size_t Widx[4];

        for(size_t kk = 0; kk<fc->nBeads; kk++)
        {
            //      for(size_t ll = kk+1; ll<fc->nBeads; ll++)
            size_t ll = kk+1; /* Note, only for elements on the first off-diagonal */
            {
                size_t idxA = kk + fc->nBeads*ll;
                if(fc->A[idxA] == 1)
                {
                    size_t N = fc->nBeads;
                    Widx[0] = kk + 2*N*ll;
                    Widx[1] = ll + 2*N*kk;

                    Widx[2] = (kk+N) + (2*N)*(ll+N);
                    Widx[3] = (ll+N) + (2*N)*(kk+N);

                    for(int ee = 0; ee<4; ee++)
                    {
                        W0[Widx[ee]] = 1;
                    }
                }
            }
        }
        return W0;
    }

    assert(0);
    return NULL;
}

int main(int argc, char ** argv)
{

    if(argc == 1)
    {
        usage();
        exit(0);
    }

    echo_args(argc, argv);

    fconf * fc = fconf_init();

    int ok = argparsing(fc, argc, argv);

    if(fc->mode == MODE_UNKNOWN)
    {
        printf("No MODE specified!\n");
        exit(1);
    }

    /* When argument == 0:
     *    Limit the memory to the amount of free memory
     *    right now
     * When argument > 0
     *     Limit memory to the number of bytes specified
     * Only implemented for linux
     */
    limit_mem(0);

    if(ok != 0)
    {
        usage();
        exit(1);
    }

    if(fc->diploid == 1) {
        printf("Diploid\n");
    } else {
        printf("Haploid\n");
    }

    /* Load the contact probability map specified with -A
     * and set fc->nBeads to size(A,1)
     */

    fconf_load_A(fc);

    /* set r0(qc) or qc(r0)
     * Can only be done when it is known how many beads
     * that are in the structures */
    set_bead_radius(fc);

    chrom * flock = load_structures(fc);

    if(fc->get_quality == 1)
    {
        if(fc->diploid == 1){ assert(0); }

        double quality = getQuality(fc, flock);
        printf("Quality: %f\n", quality);
    }

    if(fc->mode == MODE_INIT)
    {
        /* Creates initial W matrices for each structure
         * using only the contacts with theta == 1
         */
        fprintf(stdout, "Initialization mode.\n");
        fprintf(stdout, "Creating `mflock_jobs' to run with GNU parallel\n");
        FILE * jobFile = fopen("mflock_jobs", "w");
        if(jobFile == NULL)
        {
            fprintf(stderr, "Failed to create 'mflock_jobs'\n");
            exit(1);
        }

        fprintf(stdout, "Creating initial contact indication matrices, W\n");
        char * dir = malloc(32*sizeof(char));
        assert(dir != NULL);

        for(size_t kk = 0; kk<fc->nStruct; kk++)
        {
            sprintf(dir, "cf_%06zu/", kk+1);

            if(fc->rfname == NULL)
            {
                fprintf(jobFile, "mflock -w %sW.uint8.gz -o %s -x %scoords.csv --vq %f",
                        dir, dir, dir, fc->vq);
            } else {
                fprintf(jobFile, "mflock -w %sW.uint8.gz -o %s -x %scoords.csv -r %sradius.double.gz --volq %f",
                        dir, dir, dir, dir, fc->vq);
            }

            if(fc->mflock_arguments != NULL)
            {
                fprintf(jobFile, " %s", fc->mflock_arguments);
            }

            fprintf(jobFile, "\n");

            struct stat st;
            memset(&st, 0, sizeof(stat));

            if (stat(dir, &st) == -1)
            {
                mkdir(dir, 0700);
            }
        }
        free(dir);
        fclose(jobFile);

        /* Create the common start contacts */
        uint8_t * W0 = initial_W(fc);

        struct_write_W0(fc, &flock[0], W0);

        for(size_t pp = 1; pp < fc->nStruct; pp++)
        {
            printf("\r%zu", pp); fflush(stdout);
            oscp(flock[0].wfName, flock[pp].wfName);
        }
        printf("\r");

        /* Create initial radial constraints */
        if(fc->rfname != NULL)
        {
            size_t nel_r = 0;
            size_t nel_pr = 0;
            double * R = (double * ) wio_read(fc->rfname, &nel_r);
            double * PR = (double * ) wio_read(fc->prfname, &nel_pr);
            if(nel_r != nel_pr)
            {
                fprintf(stderr,
                        "ERROR: Number of elements in %s and %s does not match\n",
                        fc->rfname, fc->prfname);
                exit(-1);
            }
            if(nel_r != fc->nBeads*sizeof(double))
            {
                fprintf(stderr,
                        "Error: Number of elements in %s does not match the number of beads\n",
                        fc->rfname);
                exit(-1);
            }

            size_t used_R = 0;
            double * R0 = malloc(fc->nBeads*sizeof(double));
            assert(R0 != NULL);
            for(size_t kk = 0; kk< fc->nBeads; kk++)
            {
                R0[kk] = NAN;
                if(PR[kk] == 1)
                {
                    R0[kk] = R[kk];
                    used_R ++;
                }
            }
            printf("Using %zu elements of %s initially\n", used_R, fc->rfname);

            for(size_t pp = 0; pp < fc->nStruct; pp++)
            {
                printf("\r%zu", pp); fflush(stdout);
                sprintf(dir, "cf_%06zu/", pp+1);
                flock[pp].rfName = malloc(128*sizeof(char));
                assert(flock[pp].rfName != NULL);
                sprintf(flock[pp].rfName, "%sradius.double.gz", dir);
                struct_write_R0(fc, &flock[pp], R0);
                free(flock[pp].rfName);
                flock[pp].rfName = NULL;
            }

            printf("\r");
            free(R0);
        }

        for(size_t ff =0; ff< fc->nStruct; ff++)
            ch_free(&flock[ff]);

        free(W0);
        printf("Done!\n");
    }


    if(fc->mode == MODE_UPDATE)
    {
        fprintf(stdout, "Update/assignment mode.\n");
        /* Adds new contact to W of each structure.
         */

        double th_high = fc->th_high;
        double th_low = fc->th_low;

        if(fc->rfname != NULL)
        {
            flock_updateR(fc, flock, th_high, th_low);
        }

        /* Activation distance initialize as NAN */
        double * AD = malloc(fc->nBeads*fc->nBeads*sizeof(double));
        assert(AD != NULL);
        for(size_t kk = 0 ; kk<pow(fc->nBeads,2) ; kk++)
        {
            AD[kk] = NAN;
        }

        fprintf(stdout, "Calculating the activation distances\n");
        calc_activation_distance(fc, flock, th_high, th_low, AD);

        fprintf(stdout, "Writing to disk\n");

        pthread_t * threads = malloc(fc->nThreads*sizeof(pthread_t));
        assert(threads != NULL);
        adstruct ** S = malloc(fc->nThreads*sizeof(adstruct * ));
        assert(S != NULL);
        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            S[kk] = malloc(sizeof(adstruct));
            assert(S[kk] != NULL);
            S[kk]->thread = kk;
            S[kk]->nThreads = fc->nThreads;
            S[kk]->fc = fc;

            S[kk]->AD = AD;
            S[kk] -> th_high = th_high;
            S[kk] -> th_low = th_low;
            S[kk]->flock = flock;

            pthread_create(&threads[kk],
                           NULL,
                           run_struct_write_W_AD,
                           (void *) S[kk]);
        }

        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            pthread_join(threads[kk], NULL);
            free(S[kk]);
        }

        free(S);
        free(threads);

        printf("Done writing W-files\n");
        for(size_t kk = 0; kk<fc->nStruct; kk++)
        {
            ch_free(&flock[kk]);
        }


        printf("\r");
        printf("\n");
    }


    if(fc->mode == MODE_FINAL)
    {
        /* Loops over all the structures and
         * creates a joint contact map.
         * Note that the capture distance, cDistance is calculated
         * as a function of the bead radius. The bead radius
         * in this step does not have to match the one used in the simulations.
         */

        fprintf(stdout, ">  Finalization mode.\n");

        /* Generate contact probability matrix etc */

        const size_t nBeads = (1+fc->diploid)*fc->nBeads;
        const size_t msize = (size_t) powl(nBeads,2)*sizeof(double);
        //printf("msize=%zu\n", msize);

        pthread_t * threads = malloc(fc->nThreads*sizeof(pthread_t));
        assert(threads != NULL);
        final_tdata ** wtd = malloc(fc->nThreads*sizeof(final_tdata*));
        assert(wtd != NULL);
        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            wtd[kk] = malloc(sizeof(final_tdata));
            assert(wtd[kk] != NULL);
            wtd[kk]->thread = kk;
            wtd[kk]->nThreads = fc->nThreads;
            wtd[kk]->nBeads = (1+fc->diploid)*fc->nBeads;
            wtd[kk]->nStruct = fc->nStruct;
            wtd[kk]->M = malloc(msize);
            assert(wtd[kk]->M != NULL);
            wtd[kk]->flock = flock;
            wtd[kk]->fc = fc;

            memset(wtd[kk]->M, 0, msize);

            pthread_create(&threads[kk],
                           NULL,
                           final_tfun,
                           (void *) wtd[kk]);
        }

        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            pthread_join(threads[kk], NULL);
        }

        /* Summarize what the threads calculated from wtd */

        /* Sum up all the wanted contacts */
        double * Wtotal = malloc(pow(nBeads,2)*sizeof(double));
        assert(Wtotal != NULL);
        memset(Wtotal, 0, pow(nBeads,2)*sizeof(double));
        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            for(size_t mm = 0; mm<nBeads; mm++) {
                for(size_t nn = mm; nn<nBeads; nn++) {
                    Wtotal[mm*nBeads + nn] += (double) wtd[kk]->W[mm*nBeads + nn];
                }
            }
        }

        /* Sum up all found contacts */
        double * M = malloc(msize);
        assert(M != NULL);
        memset(M, 0, msize);
        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            for(size_t mm = 0; mm < nBeads; mm++) {
                for(size_t nn = mm; nn < nBeads; nn++) {
                    size_t idx1 = nn*nBeads + mm;
                    // size_t idx2 = mm*nBeads + nn;
                    M[idx1] += wtd[kk]->M[idx1];
                }
            }

            free(wtd[kk]->M);
        }
        /* Copy to second triangle */
        for(size_t mm = 0; mm < nBeads; mm++) {
            for(size_t nn = mm; nn < nBeads; nn++) {
                size_t idx1 = nn*nBeads + mm;
                size_t idx2 = mm*nBeads + nn;
                M[idx2] = M[idx1];
            }
        }

        /* Sum up radial profile */
        double * rprof = malloc(nBeads*sizeof(double));
        assert(rprof != NULL);
        memset(rprof, 0, nBeads*sizeof(double));
        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            for(size_t idx = 0; idx<nBeads; idx++)
            {
                rprof[idx] += wtd[kk]->rprof[idx];
            }
            free(wtd[kk]->rprof);
        }
        /* And get average */
        for(size_t idx = 0; idx<nBeads; idx++)
        {
            rprof[idx] /= fc->nStruct;
        }

        /* Free all threads and thread data remaining */
        for(size_t kk = 0; kk<fc->nThreads; kk++)
        {
            free(wtd[kk]);
        }
        free(wtd);
        free(threads);

        {
            size_t nBeads = (fc->diploid + 1)*fc->nBeads;
            for(size_t mm = 0; mm < nBeads; mm++) {
                for(size_t nn = mm; nn < nBeads; nn++) {
                    Wtotal[nn*nBeads + mm] = Wtotal[mm*nBeads + nn];
                }
            }
        }

        /* Write radial profile */
        char * rproffname = malloc(1024*sizeof(char));
        assert(rproffname != NULL);
        sprintf(rproffname, "radial_profile.csv");
        fprintf(stdout, "   Writing radial profile to: %s\n", rproffname);
        FILE * rproffile = fopen(rproffname, "w");
        for(size_t pp = 0; pp < nBeads; pp++)
        {
            fprintf(rproffile, "%f\n", rprof[pp]);
        }
        free(rprof);
        fclose(rproffile);


        /* Write contact map */
        char * mfilename = malloc(1024*sizeof(char));
        assert(mfilename != NULL);
        sprintf(mfilename, "all_contacts.double");
        fprintf(stdout, "   Writing contact map to: %s\n", mfilename);
        FILE * mout = fopen(mfilename, "w");

        if(mout == NULL)
        {
            printf("Failed to open output file %s\n", mfilename);
            exit(1);
        }

        size_t nwrite = fwrite(M, sizeof(double), (size_t) pow(nBeads,2), mout);
        printf("Wrote %zu doubles \n", nwrite);
        assert(nwrite == (size_t) pow(nBeads,2));
        fclose(mout);

        /* Write assigned contacts */
        fprintf(stdout, "   Writing sum of assigned contacts to %s\n", "Wtotal.double");
        wio_write("Wtotal.double", pow(nBeads,2)*sizeof(double), (uint8_t *) Wtotal);
        free(Wtotal);

        free(rproffname);
        free(mfilename);

        free(M);
        for(size_t kk = 0; kk<fc->nStruct; kk++)
        {
            ch_free(&flock[kk]);
        }
    }

    free(flock);
    fconf_free(fc);
    free(fc);

    return 0;
}
