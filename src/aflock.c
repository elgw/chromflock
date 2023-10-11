/**
 * @file aflock.c
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include "aflock.h"

/* @brief Limit the memory available to chromflock
 *
 * So that malloc, calloc, etc actually returns NULL at some point
 *
 * @param max_bytes
 * When argument == 0:
 *    Limit the memory to the amount of free memory
 *    right now
 * When argument > 0
 *     Limit memory to the number of bytes specified
 * Only implemented for linux
 */
static void limit_mem(size_t max_bytes)
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


static aflock * aflock_init()
{
    aflock * c = calloc(1, sizeof(aflock));
    assert(c != NULL);
    c->verbose = 1;

    c->mode = MODE_UNKNOWN;
    c->QS = 2500;

    c->ea = -1;
    c->eb = -1;
    c->ec = -1;
    c->r0 = -1;
    c->vq = -1;

    int nCpus = sysconf(_SC_NPROCESSORS_ONLN);
    /* Assuming this is the number of cores,
     * it is faster than using all threads */

    if(nCpus < 0)
    {
        printf("WARNING: Could not figure out the number of CPUs "
               "from sysconfig!\n");
        nCpus = 4;
    }

    c->nThreads = nCpus/2;
    c->nThreads == 0 ? c->nThreads = 1 : 0;

    return c;
}

static void aflock_free(aflock * c)
{
    free(c->E);
    free(c->A);
    free(c->afname);
    free(c->rfname);
    free(c->prfname);
    free(c->mflock_arguments);
}

static float eudist3(const float * A, const float * B)
{
    /* Euclidean distance between two 3D-vectors */
    return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}

static float norm3(const float * X)
{
    float n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
}

static void aflock_load_A(aflock * fc)
{

    fprintf(stdout, "Loading the Contact Probability Matrix ...\n");
    fprintf(stdout, "   '%s' ... \n", fc->afname);

    FILE * afile = fopen(fc->afname, "r");

    if(afile == NULL)
    {
        fprintf(stderr, "Failed to open file\n");
        exit(EXIT_FAILURE);
    }

    struct stat st;
    int status = stat(fc->afname, &st);
    if(status != 0)
    {
        fprintf(stderr, "Failed to get file size of %s\n", fc->afname);
        exit(EXIT_FAILURE);
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
        exit(EXIT_FAILURE);
    }

    fclose(afile);

    for(size_t kk = 0 ; kk<fsize/sizeof(double) ; kk++)
    {
        if(isfinite(fc->A[kk]))
        {
            if(fc->A[kk] < 0)
            {
                fprintf(stderr, "ERROR: A[%zu] = %f < 0\n", kk, fc->A[kk]);
                exit(EXIT_FAILURE);
            }
            if(fc->A[kk] > 1)
            {
                fprintf(stderr, "ERROR: A[%zu] = %f > 1\n", kk, fc->A[kk]);
                exit(EXIT_FAILURE);
            }
        }
    }


    fc->nBeads = sqrtl(fsize/sizeof(double));
    if(powl(fc->nBeads, 2) != fsize/sizeof(double))
    {
        fprintf(stderr, "Something is wrong with %s\n", fc->afname);
        fprintf(stderr, "The matrix does not seem to be in square form\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stdout, "   contains [%zu x %zu] elements.\n", fc->nBeads, fc->nBeads);
    return;
}

/** @brief Initialize a cf_structure
 *
 * sets the contact_pairs_file and xfName based on n
 **/
static void cf_structure_init(cf_structure * c, size_t nQ, size_t n)
{
    memset(c, 0, sizeof(cf_structure));

    c->Q = malloc(2*nQ*sizeof(uint32_t));
    assert(c->Q != NULL);

    c->contact_pairs_file = malloc(64);
    assert(c->contact_pairs_file != NULL);
    snprintf(c->contact_pairs_file, 64, "cf_%06zu/contact-pairs.u32.gz", n+1);

    c->xfName = malloc(64);
    assert(c->xfName != NULL);
    snprintf(c->xfName, 64, "cf_%06zu/coords.csv", n+1);

    return;
}

static void cf_structure_free(cf_structure * c)
{
    free(c->W);
    free(c->X);
    free(c->Q);
    free(c->contact_pairs_file);
    free(c->xfName);
    c->nQ = 0;
    return;
}


static int cf_structure_load_X(aflock * fc, cf_structure * c)
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
        exit(EXIT_FAILURE);
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
            exit(EXIT_FAILURE);
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


/** @brief Load the spatial coordinates of all structures
 *
 *
 */
static cf_structure * load_structures(aflock * fc)
{
    fprintf(stdout, "Preparing to load %zu structures ...\n", fc->nStruct);
    cf_structure * flock = malloc(fc->nStruct*sizeof(cf_structure));
    assert(flock != NULL);

    /* Load them */
    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
        printf("\r%zu", kk+1); fflush(stdout);
        cf_structure_init(&flock[kk], fc->QS, kk);
    }
    printf("\r");

    if(fc->mode == MODE_INIT)
    {
        printf("   No coordinates to be loaded\n");
    }
    else
    {
        printf("   Loading bead coordinates from %zu structures...\n",
               fc->nStruct);

        for(size_t kk = 0; kk<fc->nStruct; kk++)
        {
            printf("\r%zu", kk); fflush(stdout);
            cf_structure_load_X(fc, &flock[kk]);
        }
        printf("\r");
    }

    printf("\r");

    fprintf(stdout, "   done.\n");
    return flock;
}


static void struct_write_R0(aflock * fc, cf_structure * cc, double * R0)
{
    if(fc->diploid)
    {
        wio_write(cc->rfName, 2*fc->nBeads*sizeof(double), (void *) R0);
    } else {
        wio_write(cc->rfName, fc->nBeads*sizeof(double), (void *) R0);
    }
    return;
}


static void struct_write_W_AD(aflock * fc,
                              cf_structure * c,
                              double * AD,
                              double th_high,
                              double th_low)
{
    size_t nBeads = (1+fc->diploid)*fc->nBeads;
    size_t N = fc->nBeads;

    uint64_t nCP;
    uint32_t * CP = contact_pairs_read(c->contact_pairs_file, &nCP);
    uint8_t * W = contact_pairs_to_matrix(CP, nCP, nBeads);
    free(CP);

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

                    if(d < AD[idx1])
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

    /* Write back to disk -- TODO: only if changed ... */

    if(contact_pairs_write_from_matrix(c->contact_pairs_file, nBeads*nBeads, W))
    {
        fprintf(stderr, "Failed to write contact pairs to %s\n", c->contact_pairs_file);
        exit(EXIT_FAILURE);
    }

    free(W);
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
    printf("See the man page for more information\n");
    return;
}

static int argparsing(aflock * p, int argc, char ** argv)
{

    struct option longopts[] = {
        // MODES
        { "help",         no_argument,       NULL,  'u' },
        { "init",         no_argument,       NULL,  'I' },
        { "update",       no_argument,       NULL,  'U' },
        { "final",        no_argument,       NULL,  'F' },
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
        { "max-mem",      required_argument, NULL,  'M' },
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
                            "A:DEFIM:P:Q:R:UXh:l:n:p:r:u",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'i':
            printf("aflock (chromflock version %s)\n", cf_version);
            printf("Build date: %s, %s\n", __DATE__, __TIME__);
            printf("GIT HASH: %s\n", GIT_VERSION);
            printf("Compiler: %s\n", CC_VERSION);
            exit(EXIT_SUCCESS);
        case 'A': /* sphere contact probability matrix */
            free(p->afname);
            p->afname = malloc(strlen(optarg)+1);
            assert(p->afname != NULL);
            strcpy(p->afname, optarg);
            break;
        case 'D':
            p->diploid = true;
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
        case 'M':
            p->mem_limit = atol(optarg);
            break;
        case 'P':
            free(p->mflock_arguments);
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
            free(p->prfname);
            p->prfname = malloc(strlen(optarg)+1);
            assert(p->prfname != NULL);
            strcpy(p->prfname, optarg);
            break;
        case 'r':
            free(p->rfname);
            p->rfname = malloc(strlen(optarg)+1);
            assert(p->rfname != NULL);
            strcpy(p->rfname, optarg);
            break;
        case 'u':
            usage();
            exit(EXIT_FAILURE);
            break;
        default:
            return(EXIT_FAILURE);
        }
    }

    if(p->afname == NULL)
    {
        printf("No contact probabilty matrix (-A) file was specified\n");
        return 1;
    }

    if(p->nStruct == 0)
    {
        printf("Don't know how many structures to generate with (missing --nStruct)\n");
        return(EXIT_FAILURE);
    }

    if((p->rfname != NULL) & (p->prfname == NULL))
    {
        printf("--rpos given but not --prpos\n");
        return(EXIT_FAILURE);
    }

    if((p->rfname == NULL) & (p->prfname != NULL))
    {
        printf("--prpos given but not --rpos\n");
        return(EXIT_FAILURE);
    }

    if(p->mode == MODE_UNKNOWN)
    {
        printf("No MODE specified!\n");
        return(EXIT_FAILURE);
    }

    if(p->vq <= 0 && p->r0 <= 0)
    {
        printf("Please specify either the bead radius (--radius) or "
               "the volume quotient (--vq)\n");
        return(EXIT_FAILURE);
    }

    limit_mem(p->mem_limit);

    return EXIT_SUCCESS;
}


static void set_bead_radius(aflock * fc)
{
    double captureFactor = 4; /* TODO: deduce from dynamics program  */

    if((fc->vq > 0) && (fc->r0 > 0)) {
        fprintf(stderr, "Either -R or -Q has to be specified, not both!\n");
        exit(EXIT_FAILURE);
    }

    size_t nBeads = fc->nBeads;
    if(fc->diploid == true)
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
            exit(EXIT_FAILURE);
        }
        fc->r0 = cbrt(fc->vq*vDomain / (nBeads*4.0*M_PI/3.0));
        fc->dContact = captureFactor*fc->r0;
        allSet = 1;
    }

    if(fc->r0 > 0 && allSet == 0) {
        printf("Setting volume quotient from bead radius\n");
        if(fc->r0 > 1) {
            fprintf(stderr, "A radius of %f does not make sense\n", fc->r0);
            exit(EXIT_FAILURE);
        }
        allSet = 1;
    }

    if(allSet == 0)
    {
        fprintf(stderr, "Unable to set the bead radius. Please set either the volume quotient or the bead radius\n");
        exit(EXIT_FAILURE);
    }
    double vBeads = (double) nBeads*powl(fc->r0, 3)*4.0/3.0*M_PI;
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

/** @brief Calculate a distance threshold for each pair of beads
 *
 * This is done only for contact not already handed out. The idea is that
 * the new contacts should be handed out to the structures where the distance
 * is already the smallest under the assumption that it will be easiest for
 * those structures to integrate the new constraint.
 *
 * For high number of beads, this leads to clustering, chromosomes
 * with already many contacts get more, and those with relatively few
 * keep loosing this allocation game.
 *
*/
static void calc_activation_distance_th( adstruct * s)
{
    aflock * fc = s->fc;
    double * A = s->A;
    double * AD = s->AD;
    size_t thread = s->thread;
    size_t nThreads = s->nThreads;
    double th_low = s->th_low;
    double th_high = s->th_high;
    cf_structure * flock = s->flock;

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
                    if(fc->diploid == false) {
                        for(size_t ss = 0; ss< fc->nStruct; ss++) {
                            DS[ss] = eudist3(flock[ss].X+aa*3,
                                             flock[ss].X+bb*3);
                        }}

                    if(fc->diploid == true) {
                        for(size_t ss = 0; ss< fc->nStruct; ss++) {
                            DS[4*ss+0] = eudist3(flock[ss].X + (0+aa)*3,
                                                 flock[ss].X + (0+bb)*3);
                            DS[4*ss+1] = eudist3(flock[ss].X + (B+aa)*3,
                                                 flock[ss].X + (0+bb)*3);
                            DS[4*ss+2] = eudist3(flock[ss].X + (0+aa)*3,
                                                 flock[ss].X + (B+bb)*3);
                            DS[4*ss+3] = eudist3(flock[ss].X + (B+aa)*3,
                                                 flock[ss].X + (B+bb)*3);
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


static void * final_tfun(void * data)
{
    final_tdata * td = (final_tdata *) data;
    const size_t N = td->nBeads;

    /* This thread sums up all individual W to td->W */

    if(td->thread == 0){
        printf("   Summing up wanted contacts ... \n");
    }


    for(size_t pp = td->thread; pp < td->nStruct ; pp+=td->nThreads)
    {
        uint64_t nCP;
        uint32_t * CP = contact_pairs_read(td->flock[pp].contact_pairs_file, &nCP);
        uint8_t * W = contact_pairs_to_matrix(CP, nCP, N);
        free(CP);
        for(size_t kk = 0; kk<N*N; kk++)
        {
            td->W[kk] += W[kk];
        }
        free(W);
    }

    /* TODO:
     *
     * - dContact should be deduced from the dynamics program
     * -
     */

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

    return NULL;
}


static void * run_struct_write_W_AD(void * P)
{
    adstruct * s = (adstruct *) P;
    for(size_t kk = s->thread; kk < s->fc->nStruct ; kk += s->nThreads)
    {
        struct_write_W_AD(s->fc, &s->flock[kk], s->AD, s->th_high, s->th_low);
    }

    return NULL;
}

/* Entry point for a thread */
static void * run_calc_activation_distance_th(void * P)
{
    fflush(stdout);
    adstruct * s = (adstruct *) P;
    calc_activation_distance_th(s);
    return NULL;
}


/* Consider all positions in PR where p<th_high and p>= th_low
 * For each beads:
 * hand the constrains to the structures that has the closest radii already
 * similar to how the contact restraints are handled
 */
static void calc_activation_distance(aflock * fc, cf_structure * flock, double th_high, double th_low, double * AD)
{
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

static void flock_updateR(aflock * fc, cf_structure * flock, double th_high, double th_low)
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
            exit(EXIT_FAILURE);
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

/** @brief Hand out the first set of contacts to the structures
 *
 * Uses the contact probability matrix.
 * Any contact which has probability 1, i.e. that should appear in
 * all structures is handed out.
 */
static uint8_t * initial_W(aflock * fc)
{
    if(fc->diploid)
    {
        uint8_t * W0 = calloc(powl(2*fc->nBeads, 2), sizeof(uint8_t));
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
    } else
    {
        uint8_t * W0 = malloc(powl(fc->nBeads, 2)*sizeof(uint8_t));
        assert(W0 != NULL);
        for(size_t kk = 0; kk<powl(fc->nBeads,2); kk++)
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
}

static void chromflock_init_structures(aflock * fc, cf_structure * flock)
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
            fprintf(jobFile, "mflock --contact-pairs %scontact-pairs.u32.gz -o %s -x %scoords.csv --vq %f",
                    dir, dir, dir, fc->vq);
        } else {
            fprintf(jobFile, "mflock --contact-pairs %scontact-pairs.u32.gz -o %s -x %scoords.csv -r %sradius.double.gz --volq %f",
                    dir, dir, dir, dir, fc->vq);
        }

        if(fc->mflock_arguments != NULL)
        {
            fprintf(jobFile, " %s", fc->mflock_arguments);
        }

        fprintf(jobFile, "\n");

        struct stat st;
        memset(&st, 0, sizeof(struct stat));

        if (stat(dir, &st) == -1)
        {
            mkdir(dir, 0700);
        }
    }
    free(dir);
    fclose(jobFile);

    /* Create the common start contacts */
    uint8_t * W0 = initial_W(fc);
    uint64_t nCP = 0;
    uint32_t * CP = contact_pairs_from_matrix(W0, fc->nBeads, &nCP);
    free(W0);
    contact_pairs_write(flock[0].contact_pairs_file, CP, nCP);
    free(CP);

    for(size_t pp = 1; pp < fc->nStruct; pp++)
    {
        printf("\r%zu", pp); fflush(stdout);
        oscp(flock[0].contact_pairs_file, flock[pp].contact_pairs_file);
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
            exit(EXIT_FAILURE);
        }
        if(nel_r != fc->nBeads*sizeof(double))
        {
            fprintf(stderr,
                    "Error: Number of elements in %s does not match the number of beads\n",
                    fc->rfname);
            exit(EXIT_FAILURE);
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
        cf_structure_free(&flock[ff]);

    printf("Done!\n");
}


static void chromflock_update_structures(aflock * fc, cf_structure * flock)
{
    fprintf(stdout, "-> Update/assignment mode.\n");
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
    for(size_t kk = 0 ; kk<powl(fc->nBeads,2) ; kk++)
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
        cf_structure_free(&flock[kk]);
    }


    printf("\r");
    printf("\n");
    return;
}


/** @brief Measure the final structures
 *
 * In the 'finalization' step, some properties of the structures are measured. They are not update.
 *
 * Loops over all the structures and creates a joint contact map.
 * Note that the capture distance, cDistance is calculated as a
 * function of the bead radius. The bead radius in this step does not
 * have to match the one used in the simulations.  TODO: reduce memory
 * usage of this part
 */

static void chromflock_finalize_structures(aflock * fc, cf_structure * flock)
{
    fprintf(stdout, "Finalization mode.\n");

    /* Will not work with more threads than structures */
    fc->nThreads > fc->nStruct ? fc->nThreads = fc->nStruct : 0;

    /* Generate contact probability matrix etc */

    const size_t nBeads = (1+fc->diploid)*fc->nBeads;

    pthread_t * threads = malloc(fc->nThreads*sizeof(pthread_t));
    assert(threads != NULL);
    final_tdata ** wtd = malloc(fc->nThreads*sizeof(final_tdata*));
    assert(wtd != NULL);

    for(size_t kk = 0; kk<fc->nThreads; kk++)
    {
        wtd[kk] = calloc(1, sizeof(final_tdata));
        assert(wtd[kk] != NULL);
        wtd[kk]->thread = kk;
        wtd[kk]->nThreads = fc->nThreads;
        wtd[kk]->nBeads = nBeads;
        wtd[kk]->nStruct = fc->nStruct;
        wtd[kk]->M = calloc(nBeads*nBeads, sizeof(uint16_t));
        assert(wtd[kk]->M != NULL);
        wtd[kk]->W = calloc(nBeads*nBeads, sizeof(uint16_t));
        assert(wtd[kk]->W != NULL);
        wtd[kk]->flock = flock;
        wtd[kk]->fc = fc;

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

    /* Sum up assigned contacts */
    uint16_t * Wtotal = calloc(powl(nBeads, 2), sizeof(uint16_t));
    assert(Wtotal != NULL);
    for(size_t kk = 0; kk<fc->nThreads; kk++)
    {
        for(size_t mm = 0; mm<nBeads; mm++) {
            for(size_t nn = mm; nn<nBeads; nn++) {
                Wtotal[mm*nBeads + nn] += wtd[kk]->W[mm*nBeads + nn];
            }
        }
        free( wtd[kk]->W );
    }

    /* Copy to second triangle */
    for(size_t mm = 0; mm < nBeads; mm++) {
        for(size_t nn = mm; nn < nBeads; nn++) {
            Wtotal[nn*nBeads + mm] = Wtotal[mm*nBeads + nn];
        }
    }

    /* Write assigned contacts */
    char Wtotal_name[] = "assigned_contacts.u16";
    fprintf(stdout, "   Writing sum of assigned contacts to '%s'\n", Wtotal_name);
    FILE * fid = fopen(Wtotal_name, "w");
    if(fid == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", Wtotal_name);
        exit(EXIT_FAILURE);
    }
    size_t nWritten = fwrite(Wtotal, sizeof(uint16_t), powl(nBeads,2), fid);
    if(nWritten != powl(nBeads,2))
    {
        fprintf(stderr, "Failed to write to %s\n", Wtotal_name);
        fprintf(stderr, "Wrote %zu elements. Intended to write %zu elements\n",
                nWritten, (size_t) powl(nBeads,2));
        exit(EXIT_FAILURE);
    }
    fclose(fid);
    free(Wtotal);


    /* Sum up all found contacts based on proximity */
    uint16_t * M = calloc(nBeads*nBeads, sizeof(uint16_t));
    assert(M != NULL);
    for(size_t kk = 0; kk<fc->nThreads; kk++)
    {
        for(size_t mm = 0; mm < nBeads; mm++) {
            for(size_t nn = mm; nn < nBeads; nn++) {
                size_t idx1 = nn*nBeads + mm;
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

    /* Write contact map */
    char mfilename[] = "measured_contacts.u16";
    fprintf(stdout, "   Writing contact map to '%s'\n", mfilename);
    FILE * mout = fopen(mfilename, "w");
    if(mout == NULL)
    {
        printf("Failed to open output file %s\n", mfilename);
        exit(EXIT_FAILURE);
    }
    size_t nwrite = fwrite(M, sizeof(uint16_t), powl(nBeads,2), mout);
    if(nwrite != (size_t) powl(nBeads,2))
    {
        fprintf(stderr, "Failed to write to %s\n", mfilename);
        exit(EXIT_FAILURE);
    }
    fclose(mout);
    free(M);

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




    /* Write radial profile */
    char rproffname[] = "radial_profile.csv";
    fprintf(stdout, "   Writing radial profile to: %s\n", rproffname);
    FILE * rproffile = fopen(rproffname, "w");
    for(size_t pp = 0; pp < nBeads; pp++)
    {
        fprintf(rproffile, "%f\n", rprof[pp]);
    }
    free(rprof);
    fclose(rproffile);

    /* Free up remaining data */
    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
        cf_structure_free(&flock[kk]);
    }

    return;
}


int main(int argc, char ** argv)
{
    /* Initialize and get the default settings */
    aflock * fc = aflock_init();

    /* Parse command line arguments */
    if(argparsing(fc, argc, argv))
    {
        printf("\n");
        usage();
        exit(EXIT_FAILURE);
    }

    /* Load the contact probability map specified with -A
     * and set fc->nBeads to size(A,1)
     */

    aflock_load_A(fc);

    /* set r0(qc) or qc(r0)
     * Can only be done when it is known how many beads
     * that are in the structures */
    set_bead_radius(fc);

    /* If the structures are not loaded, this just sets the names of
     * the associates files. */
    cf_structure * flock = load_structures(fc);

    switch(fc->mode)
    {
    case MODE_INIT:
        chromflock_init_structures(fc, flock);
        break;
    case MODE_UPDATE:
        chromflock_update_structures(fc, flock);
        break;
    case MODE_FINAL:
        chromflock_finalize_structures(fc, flock);
        break;
    default:
        fprintf(stderr, "aflock does not know what to do. No mode selected\n");
    }

    free(flock);
    aflock_free(fc);
    free(fc);

    return EXIT_SUCCESS;
}
