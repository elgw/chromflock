#include "cc2cpm.h"

// TODO typedef and enum
#define MODE_DEFAULT 0
#define MODE_EQ 1

#ifndef NDEBUG
#define calloc(x,y) assert_calloc(x, y)
#define malloc(x) assert_malloc(x)
static void * assert_malloc(size_t x)
{
    double * p = (malloc)(x);
    assert(p!=NULL);
    return p;
}
static void * assert_calloc(size_t x, size_t y)
{
    double * p = (calloc)(x, y);
    assert(p!=NULL);
    return p;
}
#endif


typedef double (*statfun) (double * , size_t );

typedef struct{
    char * hFile;
    char * lFile;
    double nCont;
    size_t nStruct;
    char * aOutFile;
    char * lOutFile;
    int mode;
    int keepY;
    int verbose;
    int use_median;
} pargs;

static pargs * pargs_new();
static void pargs_show(FILE * f, pargs * p);
static void pargs_free(pargs * p);
static void toLogFile(FILE * f, int argc, char ** argv); // Write useful stuff to the log

static int argparsing(pargs * p, int argc, char ** argv);

static void help(char *);
static double double_vector_max(double * V, size_t N);
static double double_vector_mean(double * V, size_t N);



/* Remove a specific chrosmosome from H and L
 *
 * Example:
 * double * H = removeChr(H, L, &N, 23);
 * chr #23 is removed from H and L. N is updated to match the
 * new number of beads
 *
 * */
static double * removeChr(double* H0, uint8_t * L, size_t * N, uint8_t chr);

/* Write label vector */
static void
writeL(const uint8_t * L, size_t N, const char * lFile);

/*  Write contact probability matrix
 * Warning: This alters H
 */
static void writeA(double * H, size_t N2, const char * aFile);

static uint8_t * readl(const char * fileName, size_t * N); // Read label vector
static double * readh(const char * fileName, size_t * N); // Read Hi-C data

/* Get the sum of all rows (stride=1) */
//static void rsums0(double * R, // Output, vector of size N
//    double * A, // Input: Square NxN matrix
//    size_t N);
// As rsums0, however caps elements by 1
// static void rsums(double * R, double * A, size_t N);

// Max of vector of length N
static double vmax(double * V, size_t N);

/* Count the average number of contacts for a bead (defined by a vector/column, V)
 */
static double getContacts(const double * V, size_t N, double scaling, double nStruct);

/* Use binary search to figure out how to scale the Hi-C matrix to the
 * desired number of max contacts per bead
 * f determines how the scaling is evaluated, i.e., by the max contacts
 * per bead or mean ... */
static double getScaling(const double * H, const size_t N,
                         double nStruct, double nContWanted, statfun f);


/* Connect linearly consecutive beads (i.e., that are in the same chromosome) */
static void set_connectivity(double * H, uint8_t * L, size_t N, double value);

/* Get the number of neighbouring beads, i.e.
 * 0 -- only for chromosome with only one bead
 * 1 -- at the start or end of a chromosome, or
 * 2 -- inside a chromosome */
static int getNbrs(uint8_t * L, size_t N, size_t kk);

/* Figure out what bead has most contacts */
static double mostContacts(double * H, uint8_t * L, size_t N, size_t nStruct, double * cmax, size_t * nmax);

/* Some sanity checking of the Hi-C matrix */
static void checkH(double * H, size_t N);

/* Set diagonal to 0
 * diag=0 -- only the main diagonal
 * diag=1 -- also the +/- 1 diagonals, etc */
static void clearDiagonal(double * H, size_t N, int diag);
static double vmax(double * V, size_t N);

/* Clears beads for which look suspect in the Hi-C matrix having more contacts to some other bead than to itself */
static void clearWeak(double * H, size_t N);

/* Scale H to have at most nCont contacts per bead */
static int scale(double * H, uint8_t * L, size_t N,
                 size_t nStruct, double nCont, statfun stat);


static pargs * pargs_new()
{
    pargs * args = calloc(1, sizeof(pargs));
    /* HiC data, encodes as double */
    args->hFile = NULL;
    /* Labels, encoded as uint8 */
    args->lFile = NULL;
    /* Contacts per bead */
    args->nCont = 8;
    /* Number of structures to be generated */
    args->nStruct = 10000;
    /* Contact probability matrix, output file name */
    args->aOutFile = NULL;
    args->lOutFile = NULL;
    /* Normalization or not? */
    args->mode = MODE_DEFAULT;
    /* Remove chrY from the output matrix? */
    args->keepY = 0;
    args->verbose = 1;
    args->use_median = 0;
    return args;
}

static int argparsing(pargs * p, int argc, char ** argv)
{

    struct option longopts[] = {
        {"aOut",    required_argument, NULL, 'A'},
        {"nCont",   required_argument, NULL, 'c'},
        {"mode_eq", no_argument,       NULL, 'e'},
        {"help",    no_argument,       NULL, 'h'},
        {"mean",    no_argument,       NULL, 'm'},
        {"hFile",   required_argument, NULL, 'H'},
        {"lFile",   required_argument, NULL, 'L'},
        {"nStruct", required_argument, NULL, 's'},
        {"version", no_argument,       NULL, 'v'},
        {"y",       no_argument,       NULL, 'y'},
        {NULL,      0,                 NULL,  0 }
    };

    int ch;
    while( ( ch = getopt_long(argc, argv, "A:c:ehH:L:ms:vy", longopts, NULL)) != -1)
    {
        switch (ch)
        {
        case 'A':
            p->aOutFile = strdup(optarg);
            break;
        case 'c':
            p->nCont = atof(optarg);
            break;
        case 's':
            p->nStruct = atol(optarg);
            break;
        case 'e':
            p->mode = MODE_EQ;
            break;
        case 'm':
            p->use_median = 1;
            break;
        case 'y':
            p->keepY = 1;
            break;
        case 'L':
            p->lFile = realpath(optarg, NULL);
            if(p->lFile == NULL)
            {
                printf("Can't find the Label file\n");
                exit(1);
            }
            break;
        case 'H':
            p->hFile = realpath(optarg, NULL);
            if(p->hFile == NULL)
            {
                printf("Can't find the Hi-C file\n");
                exit(1);
            }
            break;
        case 'h':
            help(argv[0]);
            exit(0);
        case 'v':
            printf("cc2cpm (chromflock version %s)\n", cf_version);
            printf("Build date: %s, %s\n", __DATE__, __TIME__);
            printf("GIT HASH: %s\n", GIT_VERSION);
            printf("Compiler: %s\n", CC_VERSION);
            exit(0);
        default:
            break;
        }
    }
    int ret = 0;
    if(p->lFile == NULL)
    {
        printf("--lFile not specified\n");
        ret = 1;
    }
    if(p->hFile == NULL)
    {
        printf("--hFile not specified\n");
        ret = 1;
    }

    if(ret != 0)
        return ret;

    char ystr[] = "_with_ChrY";
    if(p->keepY != 1)
    {
        sprintf(ystr, "%s", "");
    }
    char * modeStr = malloc(100*sizeof(char));
    sprintf(modeStr, "%s", "");


    if(p->mode == MODE_EQ)
    {
        if(p->use_median)
        {
            sprintf(modeStr, "%s", "_EQ_MEAN");
        } else {
            sprintf(modeStr, "%s", "_EQ_MAX");
        }
    }
    if(p->mode == MODE_DEFAULT)
    {
        if(p->use_median)
        {
            sprintf(modeStr, "%s", "_MEAN");
        } else {
            sprintf(modeStr, "%s", "_MAX");
        }
    }

    if(p->aOutFile == NULL)
    {
        p->aOutFile = malloc(1024*sizeof(char));
        sprintf(p->aOutFile, "%f_%zu%s%s.A.double", p->nCont, p->nStruct, modeStr, ystr);
        p->lOutFile = malloc(1024*sizeof(char));
        sprintf(p->lOutFile, "%f_%zu%s%s.L.uint8", p->nCont, p->nStruct, modeStr, ystr);
    } else {
        p->lOutFile = malloc(1024*sizeof(char));
        sprintf(p->lOutFile, "%s.L.uint8", p->aOutFile);
    }
    return ret;
}

static double min_double(double a, double b)
{
    if(a<b)
    {
        return a;
    }
    return b;
}

static void writeL(const uint8_t * L, size_t N, const char * lFile)
{
    FILE * f = fopen(lFile, "wb");
    if(f == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", lFile);
        exit(EXIT_FAILURE);
    }
    size_t nwritten = fwrite(L, sizeof(uint8_t), N, f);
    fclose(f);
    if(nwritten != N)
    {
        fprintf(stderr, "Could not write %s to disk\n", lFile);
        fprintf(stderr, "Tried to write %zu, but could only write %zu\n",
                N, nwritten);
        exit(EXIT_FAILURE);
    }
    return;
}

static void writeA(double * H, size_t N2, const char * aFile)
{

    for(size_t kk = 0; kk<N2; kk++)
    {
        H[kk] = min_double(H[kk], 1);
    }

    FILE * f = fopen(aFile, "wb");
    if(f == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", aFile);
        exit(EXIT_FAILURE);
    }
    size_t nwritten = fwrite(H, sizeof(double), N2, f);
    if(nwritten != N2)
    {
        fprintf(stderr, "Could not write %s to disk\n", aFile);
        fprintf(stderr, "Tried to write %zu, but could only write %zu\n",
                N2, nwritten);
        exit(EXIT_FAILURE);
    }
    fclose(f);
    return;
}

#if 0
static void rsums(double * R, double * A, size_t N)
{
    // Sum of all rows
    // However does not count any element as more than 1.
    //
    for(size_t kk = 0; kk<N; kk++)
    {
        R[kk] = 0;
        for(size_t ll = 0; ll<N; ll++)
        {
            double v = min_double(A[kk*N + ll], 1);
            R[kk] += v;
        }
    }
}
#endif

/* Get number of contacts for a bead (average for the nStruct)
 * specified by a row from the HiC-matrix, V, with N elements
 * More specific, scaling*V is used
 */

static double getContacts(const double * V, size_t N, double scaling, double nStruct)
{
    double v = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        v += round(scaling*V[kk]*nStruct);
    }
    return v/nStruct;
}

static double double_vector_max(double * V, size_t N)
{
    double max = V[0];
    for(size_t kk = 0; kk<N; kk++)
    {
        V[kk] > max ? max = V[kk] : 0;
    }
    return max;
}


static double double_vector_mean(double * V, size_t N)
{
    double sum = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        sum += V[kk];
    }
    return sum / (double) N;
}

static double nContacts(const double * H, const size_t N, double scaling, double nStruct, statfun f)
{
    double * nCont = malloc(N*sizeof(double));
    for(size_t kk = 0; kk<N; kk++)
    {
        nCont[kk] = getContacts(H+kk*N, N, scaling, nStruct);
    }
    double v = f(nCont, N);
    free(nCont);
    return v;
}

static double getScaling(const double * H, const size_t N,
                         double nStruct, double nContWanted,
                         statfun stat)
{
    double args[3] = {0, 50, 100};
    double vals[3];

    vals[0] = nContacts(H, N, args[0], nStruct, stat);
    vals[2] = nContacts(H, N, args[2], nStruct, stat);

    while(vals[2] < nContWanted)
    {
        args[2] *= 2.0;
        vals[2] = nContacts(H, N, args[2], nStruct, stat);
    }

    args[1] = 0.5*(args[0] + args[2]);
    vals[1] = nContacts(H, N, args[1], nStruct, stat);
#if 0
    printf("Starting point:\n");
    printf("args = [%f, %f, %f]\n", args[0], args[1], args[2]);
    printf("vals = [%f, %f,  %f]\n", vals[0], vals[1], vals[2]);
#endif

    size_t iter = 1;
    while( (fabs(vals[1] - nContWanted) > 1e-3)
           && (args[2]-args[0] > 1e-6))
    {
        if(vals[1] > nContWanted)
        {
            vals[2] = vals[1];
            args[2] = args[1];
        } else {
            vals[0] = vals[1];
            args[0] = args[1];
        }
        args[1] = 0.5*(args[0]+args[2]);
        vals[1] = nContacts(H, N, args[1], nStruct, stat);
#if 0
        printf("iter: %zu\n", iter);
        printf("args = [%f, %f, %f]\n", args[0], args[1], args[2]);
        printf("vals = [%f, %f,  %f]\n", vals[0], vals[1], vals[2]);
        printf("\n");
#endif
        iter++;
    }
#if 0
    printf("Final contacts per bead: %f (wanted %f) for scaling=%f\n", vals[1], nContWanted, args[1]);
    printf("args = [%f, %f, %f]\n", args[0], args[1], args[2]);
    printf("vals = [%f, %f,  %f]\n", vals[0], vals[1], vals[2]);
#endif
    printf("Found the scaling factor: %f in %zu iterations\n", args[1], iter);
    return args[1];
}


static void set_connectivity(double * H, uint8_t * L, size_t N, double value)
/* Set the matrix elements connecting consecutive beads in H to 0
 * Other elements are left untouched
 * They will later on be counted as 1 regardless of scale
 */
{
    for(size_t kk = 1; kk < N; kk++)
    {
        if(L[kk] == L[kk-1])
        {
            size_t a = kk;
            size_t b = kk-1;
            H[N*a + b] = value;
            H[N*b + a] = value;
        }
    }
}


static uint8_t * readl(const char * fileName, size_t * N)
{
    fprintf(stdout, " >> Loading  %s  \n", fileName);

    FILE * hfile = fopen(fileName, "r");

    if(hfile == NULL)
    {
        fprintf(stderr, "    ! Failed to open %s\n", fileName);
        exit(EXIT_FAILURE);
    }

    struct stat st;
    int status = stat(fileName, &st);
    if(status != 0)
    {
        fprintf(stderr, "    ! Failed to get file size of %s\n", fileName);
        exit(EXIT_FAILURE);
    }

    size_t fsize = st.st_size;

    uint8_t * L = malloc(fsize*sizeof(uint8_t));
    if(L == NULL)
    {
        fprintf(stderr, "    ! Failed to allocate memory for L\n");
        exit(-1);
    }

    size_t nRead = fread(L, sizeof(uint8_t), fsize/sizeof(uint8_t), hfile);

    if(nRead*sizeof(uint8_t) != fsize)
    {
        fprintf(stderr, "    ! nRead = %zu, fsize = %zu\n", nRead, fsize);
        exit(-1);
    }

    N[0] = fsize/sizeof(uint8_t);
    fclose(hfile);
    return L;

}

static double * readh(const char * fileName, size_t * _N)
{
    fprintf(stdout, " >> Loading  %s  \n", fileName);

    FILE * hfile = fopen(fileName, "r");

    if(hfile == NULL)
    {
        fprintf(stderr, "    ! Failed to open file\n");
        exit(EXIT_FAILURE);
    }

    struct stat st;
    int status = stat(fileName, &st);
    if(status != 0)
    {
        fprintf(stderr, "    ! Failed to get file size of %s\n", fileName);
        exit(EXIT_FAILURE);
    }

    size_t fsize = st.st_size;

    double * A = malloc(fsize);
    if(A == NULL)
    {
        fprintf(stderr, "    ! Failed to allocate memory for A\n");
        exit(EXIT_FAILURE);
    }

    size_t nRead = fread(A, sizeof(double), fsize/sizeof(double), hfile);

    if(nRead*sizeof(double) != fsize)
    {
        fprintf(stderr, "    ! nRead = %zu, fsize = %zu\n", nRead, fsize);
        exit(EXIT_FAILURE);
    }

    size_t N = fsize/sizeof(double);
    fclose(hfile);

    // Convert nan (and, any inf present?) to 0
    int onlyNormal = 1;
    for(size_t kk = 0; kk<N; kk++)
    {
        if(isnormal(A[kk]) != 1)
        {
            A[kk] = 0;
            onlyNormal = 0;
        }
    }
    if(onlyNormal == 0)
    {
        printf("Warning: Non-finite values encountered. "
               "Those were converted to 0\n");
    }

    size_t SN = sqrt(N);
    if(SN*SN != N)
    {
        fprintf(stderr, "%s does not seem to be a square matrix\n", fileName);
        fprintf(stderr, "it has %zu elements\n", N);
        exit(EXIT_FAILURE);
    }

    _N[0] = SN;

    return A;

}


static int getNbrs(uint8_t * L, size_t N, size_t kk)
{
    /* Add the contacts to possible neighbours in the same chromosome
     */
    int c = 0;
    if(kk>0)
    {
        if(L[kk-1] == L[kk])
        {
            c++;
        }}
    if(kk+1<N)
    {
        if(L[kk] == L[kk+1])
        {
            c++;
        }}
    return c;
}

static void help(char * myname)
{
    printf("Usage:\n");
    printf("%s --hFile H.double --lFile L.uint8 --nCont n --nStruct n\n", myname);
    printf("Required arguments:\n");
    printf(" --hFile <file>\n\t HiC-matrix encoded as raw double\n");
    printf(" --lFile <file>\n\t Chromosome label vector encoded as raw uint8_t\n");
    printf("Optional arguments:\n");
    printf(" --nStruct n\n\t Number of structures that you will ask chromflock "
           "to generate later on\n");
    printf(" --nCont n\n\t max number of contacts per bead (for the specific "
           "number of structures)\n");
    printf(" --mean\n\t"
           "scale the matrix so that the mean number of contacts per bead is nCont\n");
    printf(" --aOut <file>\n\t name of output contact probability matrix\n");
    printf(" --mode_eq\n\t Give equal number of contacts for each bead by "
           "balancing without the -1, 0, 1 diagonals. The balancing does not "
           "care about rounding so only in the limit case when nStruct is large "
           "will a constant number of contacts per bead be reached.\n");
    printf(" --y\n\t Keep chrY encoded by 24 (by default it is removed)\n");
    printf(" --help\n\t Show a short help section\n");
    printf("\n");
    printf("see the man page for more information\n");
    exit(0);
}

static double mostContacts(double * H, uint8_t * L, size_t N, size_t nStruct, double * cmax, size_t * nmax)
{
    /* Count number of contacts per bead if L is supplied the
     * connection between adjacent beads is assumed to be set to 0 and
     * counts are increased based on L */

    cmax[0] = 0;
    nmax[0] = 0;
    double ctot = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        double nContacts = getContacts(H+kk*N, N, 1, nStruct);
        if(L != NULL) {
            nContacts += getNbrs(L, N, kk);
        }

        ctot += nContacts;
        if(nContacts > cmax[0])
        {
            cmax[0] = nContacts;
            nmax[0] = kk;
        }
    }
    return ctot/N;
}

static void checkH(double * H, size_t N)
{
    for(int kk = 0 ; kk<10 ; kk++)
    {
        printf("%e ", H[kk]);
    }
    printf("\n");

    for(size_t kk = 0; kk<N*N; kk++)
    {
        if(isnormal(H[kk]) != 1)
        {
            H[kk] = 0;
        }
    }

    double hmin = 1e99; double hmax = 0;
    for(size_t kk = 0; kk<N*N; kk++)
    {
        if(H[kk] > hmax)
            hmax = H[kk];
        if(H[kk] < hmin)
            hmin = H[kk];
    }
    printf("   H is in [%f, %f]\n", hmin, hmax);

    int isSymmetric = 1;
    for(size_t nn = 0; nn<N; nn++)
    {
        for(size_t mm = 0; mm<N; mm++)
        {
            if( fabs(H[nn + N*mm] - H[mm + N*nn]) > 1e-9 )
            {
                isSymmetric = 0;
            }
        }
    }
    if(isSymmetric == 1){
        printf("    H^T == H\n");
    } else {
        printf(" ERROR: H^T != H\n");
        printf("        Please check the input data!\n");
        printf("        Press Ctrl+C to abort (or enter to continue if you dare)\n");
        getchar();
    }

    double diagsum = 0;
    for(size_t kk = 0; kk<N*N; kk+= N+1)
    { diagsum+=H[kk]; }

    if(diagsum == 0)
    {
        printf("ERROR: Sum of diagonal elements is 0! This can't raw Hi-C or TCC data.\n");
        exit(1);
    }


    return;
}

static void clearDiagonal(double * H, size_t N, int diag)
{
    /* Set diagonal(s) to 0 , diag=0 corresponds to the main diagonal
     * diag=1 refers to the +/- 1 diags, etc, ... */
    for(int dd = 0 ; dd<=diag; dd++)
    {
        for(size_t kk = dd; kk<N*N; kk+=(N+1))
        {
            H[kk] = 0;
        }
        if(dd != 0)
        {
            for(size_t kk = N*dd; kk<N*N; kk+=(N+1))
            {
                H[kk] = 0;
            }
        }
    }
}

static double vmax(double * V, size_t N)
/* Return largest value of a vector. Only normal values are considered
 */
{
    double ma = 1e-99;
    for(size_t kk = 0; kk<N; kk++)
    {
        if(isnormal(V[kk]))
        {
            if(V[kk] > ma)
            {
                ma = V[kk];
            }
        }
    }
    return ma;
}

static void clearWeak(double * H, size_t N)
{
    for(size_t dd = 0; dd<N; dd++)
    {
        if(H[dd + N*dd] < vmax(H+dd*N, N))
        {
            for(size_t mm = 0; mm<N; mm++)
                H[mm + N*dd] = 0;
            for(size_t mm = 0; mm<N; mm++)
                H[dd + N*mm] = 0;
        }
    }
    return;
}

static void pargs_show(FILE * f, pargs * p)
{
    fprintf(f, " Hi-C matrix:\n\t %s\n", p->hFile);
    fprintf(f, " Chromosome labels:\n\t %s\n", p->lFile);
    fprintf(f, " Number of structures: %zu\n", p->nStruct);
    if(p->use_median) {
        fprintf(f, " Wanted MEAN contacts per bead: %f\n", p->nCont);
    } else {
        fprintf(f, " Wanted MAX contacts per bead: %f\n", p->nCont);
    }
    fprintf(f, " OUT: Contact probability matrix:\n\t %s\n", p->aOutFile);
    fprintf(f, " OUT: Chromosome label vector:\n\t %s\n", p->lOutFile);
    fprintf(f, " MODE: ");
    switch(p->mode)
    {
    case MODE_DEFAULT:
        fprintf(f, "DEFAULT -- scale to max wanted contacts, set connectivity\n");
        break;
    case MODE_EQ:
        fprintf(f, "EQ -- equal number of contacts per bead\n");
        break;
    }

    if(p->keepY == 0)
    {
        fprintf(f, " chrY will be removed\n");
    } else {
        fprintf(f, " chrY will NOT be removed\n");
    }

    return;
}

static void pargs_free(pargs * p)
{
    // Free up ...
    if(p->hFile != NULL)
        free(p->hFile);
    if(p->lFile != NULL)
        free(p->lFile);
    if(p->aOutFile != NULL)
        free(p->aOutFile);
    if(p->lOutFile != NULL)
        free(p->lOutFile);
    free(p);
}

static void toLogFile(FILE * f, int argc, char ** argv)
{
    fprintf(f, "CC_VERSION: %s\n", CC_VERSION);
    fprintf(f, "GIT_VERSION: %s\n", GIT_VERSION);

    fprintf(f, "Command line:\n");
    for(int kk = 0; kk<argc; kk++)
    {
        fprintf(f, "%s ", argv[kk]);
    }
    fprintf(f, "\n");

    time_t current_time;
    char* c_time_string;

    /* Obtain current time. */
    current_time = time(NULL);

    if (current_time == ((time_t)-1))
    {
        (void) fprintf(stderr, "Failure to obtain the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    if (c_time_string == NULL)
    {
        (void) fprintf(stderr, "Failure to convert the current time.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "Current time is %s", c_time_string);
}

static double * removeChr(double* H0, uint8_t * L, size_t * N, uint8_t chr)
{

    size_t N0 = N[0];
    // Count the number of beads with the chr to be excluded
    size_t nx = 0;
    for(size_t kk = 0; kk<N0; kk++)
    {
        if(L[kk] == chr)
        {
            nx++;
        }
    }

    if(nx == 0)
    {
        //printf("Nothing to do\n");
        return H0;
    }

    size_t N1 = N0-nx;

    // Create new H
    double * H = malloc(pow(N0-nx, 2)*sizeof(double));
    size_t nnw = 0;
    for(size_t nn = 0; nn < N0; nn++)
    {
        if(L[nn] == chr)
        { nnw++; }
        else
        {
            size_t mmw = 0;
            for(size_t mm = 0; mm < N0; mm++)
            {
                if(L[mm] == chr)
                { mmw++; }
                else
                {
                    H[N1*(nn-nnw) + (mm-mmw)] = H0[N0*nn + mm];
                }
            }
        }
    }
    free(H0);

    // Update L
    size_t writepos = 0;
    for(size_t kk = 0; kk<N0; kk++)
    {
        if(L[kk] != chr)
        {
            L[writepos++] = L[kk];
        }
    }

    N[0] = N1;
    return H;
}

static int scale(double * H, uint8_t * L, size_t N, size_t nStruct, double nCont, statfun stat)
{
    /* Remove connection between consecutive beads those values can't be scaled */
    set_connectivity(H, L, N, 0.0);

    double scale = getScaling(H, N, nStruct, nCont-2.0, stat);
    for(size_t kk = 0; kk<N*N; kk++)
    {
        H[kk] *= scale;
    }

    // Put connectivity back
    set_connectivity(H, L, N, 1.0);

    return 0;
}

int cc2cpm(int argc, char ** argv)
{
    pargs * args = pargs_new();
    if( argparsing(args, argc, argv) )
    {
        help(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Set number of contacts to mean or max over the beads? */
    statfun stat = double_vector_max;
    if(args->use_median)
    {
        stat = double_vector_mean;
    }

    pargs_show(stdout, args);
    printf("Press <Enter> to continue\n");
    getchar();

    size_t nL = 0;
    uint8_t * L = readl(args->lFile, &nL);
    uint8_t lmin = 255; uint8_t lmax = 0;
    for(size_t kk = 0; kk<nL; kk++)
    {
        if(lmin > L[kk])
            lmin = L[kk];
        if(lmax < L[kk])
            lmax = L[kk];
    }
    printf("    L has %zu elements, min: %u, max %u\n", nL, lmin, lmax);

    size_t N = 0;
    double * H = readh(args->hFile, &N);
    printf("    H has %zu elements, [%zu x %zu]\n", N*N, N, N);

    if(N != nL)
    {
        fprintf(stderr, "The number of elements in the Labels file does not match the Hi-C matrix size\n");
        fprintf(stderr, "%zu vs %zu\n", nL, N);
        exit(EXIT_FAILURE);
    }

    size_t nmax = 0;
    double cmax = -1e99;
    double cmean = mostContacts(H, NULL, N, args->nStruct, &cmax, &nmax);

    printf(" >> Before any processing:\n");
    printf("    Most contacts for bead/bin %zu : %f\n", nmax, cmax);
    printf("    Mean number of contacts: %f\n", cmean);
    checkH(H, N);

    if(args->keepY == 0)
    {
        printf(" >> Removing ChrY\n");
        size_t N0 = N;
        H = removeChr(H, L, &N, 24);
        if(N == N0)
        {
            printf("    There was no ChrY to remove!\n");
        }
    }

    printf(" >> Removing contacts for beads where N(i,i) != max(N(i,:))\n");
    clearWeak(H, N);

    printf("    Setting diag(H) = 0\n");
    clearDiagonal(H, N, 0);

    if(args->mode == MODE_DEFAULT)
    {
        printf(" >> Finding scaling factor\n");
        scale(H, L, N, args->nStruct, args->nCont, stat);
    }

    if(args->mode == MODE_EQ)
    {
        printf(" >> Balancing H\n");
        double berror = balance(H, N);
        printf("    Largest error: %f\n", berror);
        //  writeA(H, N2, "H_balanced.double");

        printf(" >> Finding scaling factor\n");
        scale(H, L, N, args->nStruct, args->nCont, stat);
    }

    printf(" >> Writing ...\n");
    printf("    Press enter to continue, or abort with Ctrl+C\n");
    getchar();

    char * logFileName = malloc(strlen(args->aOutFile) + 1024);
    sprintf(logFileName, "cc2cpm.log");


    if(1){
        FILE * logFile = fopen(logFileName, "w");
        if(logFile == NULL)
        {
            printf("Failed to open log file : %s\n", logFileName);
        } else {
            toLogFile(logFile, argc, argv);
            pargs_show(logFile, args);
            fclose(logFile);
        }
    }

    writeA(H, N*N, args->aOutFile);
    writeL(L, N,  args->lOutFile);

    free(H);
    free(L);
    pargs_free(args);
    free(logFileName);
    printf(" >> Done\n");
    return 0;
}
