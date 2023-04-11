#include "cc2cpm.h"

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
} pargs;

static pargs * pargs_new();
static void pargs_show(FILE * f, pargs * p);
static void pargs_free(pargs * p);
static void toLogFile(FILE * f, int argc, char ** argv); // Write useful stuff to the log

static int argparsing(pargs * p, int argc, char ** argv);

static void help(char *);


/* Remove a specific chrosmosome from H and L
 *
 * Example:
 * double * H = removeChr(H, L, &N, 23);
 * chr #23 is removed from H and L. N is updated to match the
 * new number of beads
 *
 * */
static double * removeChr(double* H0, uint8_t * L, size_t * N, uint8_t chr);

/* Read and write as uint8_t and double */
static void writeL(uint8_t * L, size_t N, char * lFile); // Write label vector
static void writeA(double * H, size_t N2, char * aFile); // Write contact probability matrix
static uint8_t * readl(char * fileName, size_t * N); // Read label vector
static double * readh(char * fileName, size_t * N); // Read Hi-C data

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
static double getContacts(double * V, size_t N, double scaling, double nStruct);

/* Use binary search to figure out how to scale a column to the desired number of
 * contacts*/
static double getScaling(double * V, size_t N, double nCont, double nStruct, double maxScaling);

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
static int scale(double * H, uint8_t * L, size_t N, size_t nStruct, double nCont);


static pargs * pargs_new()
{
    pargs * args = malloc(sizeof(pargs));
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
    return args;
}

static int argparsing(pargs * p, int argc, char ** argv)
{

    struct option longopts[] = {
        {"aOut",    required_argument, NULL, 'A'},
        {"nCont",   required_argument, NULL, 'c'},
        {"mode_eq", no_argument,       NULL, 'e'},
        {"help",    no_argument,       NULL, 'h'},
        {"hFile",   required_argument, NULL, 'H'},
        {"lFile",   required_argument, NULL, 'L'},
        {"nStruct", required_argument, NULL, 's'},
        {"version", no_argument,       NULL, 'v'},
        {"y",       no_argument,       NULL, 'y'},
        {NULL,      0,                 NULL,  0 }
    };

    int ch;
    while( ( ch = getopt_long(argc, argv, "A:c:ehH:L:s:vy", longopts, NULL)) != -1)
    {
        switch (ch)
        {
        case 'c':
            p->nCont = atof(optarg);
            break;
        case 's':
            p->nStruct = atol(optarg);
            break;
        case 'e':
            p->mode = MODE_EQ;
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
        sprintf(modeStr, "%s", "_EQ");
    }
    if(p->mode == MODE_DEFAULT)
    {
        sprintf(modeStr, "%s", "_MAX");
    }

    if(p->aOutFile == NULL)
    {
        p->aOutFile = malloc(1024*sizeof(char));
        sprintf(p->aOutFile, "%f_%zu%s%s.A.double", p->nCont, p->nStruct, modeStr, ystr);
        p->lOutFile = malloc(1024*sizeof(char));
        sprintf(p->lOutFile, "%f_%zu%s%s.L.uint8", p->nCont, p->nStruct, modeStr, ystr);
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

#if 0
static double max_double(double a, double b)
{
    if(a>b)
    {
        return a;
    }
    return b;
}
#endif

static void writeL(uint8_t * L, size_t N, char * lFile)
{
    FILE * f = fopen(lFile, "wb");
    fwrite(L, sizeof(uint8_t), N, f);
    fclose(f);
    return;
}

static void writeA(double * H, size_t N2, char * aFile)
{

    for(size_t kk = 0; kk<N2; kk++)
    {
        H[kk] = min_double(H[kk], 1);
    }

    FILE * f = fopen(aFile, "wb");
    fwrite(H, sizeof(double), N2, f);
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

static double getContacts(double * V, size_t N, double scaling, double nStruct)
{
    // Get number of contacts for a bead (average for the nStruct)
    // specified by a row from the HiC-matrix, V, with N elements
    // More specific, scaling*V is used

    double v = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        v += round(scaling*V[kk]*nStruct);
    }
    return v/nStruct;
}

static double getScaling(double * V, size_t N, double nCont, double nStruct, double maxScaling)
/* Using binary search to find the optimal scaling factor */
{
    double args[3] = {1e-9, 50, 100};
    args[2] = maxScaling;
    args[1] = (args[0]+args[2])/2.0;
    double vals[3];

    for(int kk = 0; kk<3; kk++)
    {
        vals[kk] = getContacts(V, N, args[kk], nStruct);
        //    printf("vals[%d] = %f\n", kk, vals[kk]);
    }

    assert(vals[0]-nCont < 0);
    assert(vals[2]-nCont > 0);

    // While search interval is large enough
    while((args[2] - args[0]) > 1e-11)
    {
        if((vals[1]-nCont) < 0)
        {
            args[0] = args[1];
            vals[0] = vals[1];
        } else {
            args[2] = args[1];
            vals[2] = vals[1];
        }

        args[1] = (args[0] + args[2])/2.0;
        vals[1] = getContacts(V, N, args[1], nStruct);

        if(0){
            for(int kk = 0; kk<3; kk++)
            {
                printf("args[%d] = %f vals[%d] = %f\n", kk, args[kk], kk, vals[kk]);
            }
            printf("\n");
        }
    }
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


static uint8_t * readl(char * fileName, size_t * N)
{
    fprintf(stdout, " >> Loading  %s  \n", fileName);

    FILE * hfile = fopen(fileName, "r");

    if(hfile == NULL)
    {
        fprintf(stderr, "    ! Failed to open file\n");
        exit(-1);
    }

    struct stat st;
    int status = stat(fileName, &st);
    if(status != 0)
    {
        fprintf(stderr, "    ! Failed to get file size of %s\n", fileName);
        exit(-1);
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

static double * readh(char * fileName, size_t * N)
{
    fprintf(stdout, " >> Loading  %s  \n", fileName);

    FILE * hfile = fopen(fileName, "r");

    if(hfile == NULL)
    {
        fprintf(stderr, "    ! Failed to open file\n");
        exit(-1);
    }

    struct stat st;
    int status = stat(fileName, &st);
    if(status != 0)
    {
        fprintf(stderr, "    ! Failed to get file size of %s\n", fileName);
        exit(-1);
    }

    size_t fsize = st.st_size;

    double * A = malloc(fsize);
    if(A == NULL)
    {
        fprintf(stderr, "    ! Failed to allocate memory for A\n");
        exit(-1);
    }

    size_t nRead = fread(A, sizeof(double), fsize/sizeof(double), hfile);

    if(nRead*sizeof(double) != fsize)
    {
        fprintf(stderr, "    ! nRead = %zu, fsize = %zu\n", nRead, fsize);
        exit(-1);
    }

    N[0] = fsize/sizeof(double);
    fclose(hfile);

    // Convert nan (and, any inf present?) to 0
    int onlyNormal = 1;
    for(size_t kk = 0; kk<N[0]; kk++)
    {
        if(isnormal(A[kk]) != 1)
        {
            A[kk] = 0;
            onlyNormal = 0;
        }
    }
    if(onlyNormal == 0)
    {
        printf("Warning: Non-finite values encountered. Those were converted to 0\n");
    }
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
    /* Count number of contacts per bead
     * if L is supplied the connection between adjacent beads is assumed to be set to 0
     * and counts are increased based on L */
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
    fprintf(f, " Wanted max contacts per bead: %f\n", p->nCont);
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

static int scale(double * H, uint8_t * L, size_t N, size_t nStruct, double nCont)
{
    size_t iter = 0;
    size_t iter_max = 100;

    double cmax = 0;
    size_t nmax = 0;

    /* Remove connection between consecutive beads those values can't be scaled */
    set_connectivity(H, L, N, 0.0);

    /* Find the bead with most contacts */
    mostContacts(H, L, N, nStruct, &cmax, &nmax);
    printf("    Using bead %zu with %f contacts/%f for initial attempt\n",
           nmax, cmax, nCont);

    double max_scale = 1000;
    do
    {
        // Get optimal scaling factor based on a specific bead
        // it would be better to do it on all at the same time but it
        // is of course faster to do it on a single one
        //
        double scale = getScaling(H+N*nmax, N, nCont - getNbrs(L, N, nmax), nStruct, max_scale);
        max_scale = 1+1e-4; //
        printf("    Scaling factor: %f\n", scale);


        // Use the found scaling factor
        for(size_t kk = 0; kk<N*N; kk++)
        {
            H[kk] *= scale;
        }
        //   double nGot = getContacts(H+N*nmax, N, 1, nStruct);
        //   double nBrs = getNbrs(L, N, nmax);
        //    printf(" ---> nGot = %f+%f, N = %zu\n", nGot, nBrs, N);


        // And see if it works out fine
        // For low number of structures the rounding procedure sometimes
        // cause the initial guess to be wrong
        // No guarantee that it will converge
        printf("     >> Validating\n");

        double ctot = mostContacts(H, L, N, nStruct, &cmax, &nmax);

        printf("        Found most contacts for bead %zu: %f\n", nmax, cmax);
        printf("        Average: %f contacts per bead\n", ctot);
        iter++;
    } while( fabs(cmax - nCont) > 0.01 &&
             iter < iter_max );

    if(iter == iter_max)
    {
        printf("Max number of iterations reached\n");
        exit(1);
    }

    // Put connectivity back
    set_connectivity(H, L, N, 1.0);

    return 0;
}

int cc2cpm(int argc, char ** argv)
{

    pargs * args = pargs_new();
    if( argparsing(args, argc, argv))
    {
        help(argv[0]);
        exit(-1);
    }

    char * hFile = args->hFile;
    char * lFile = args->lFile;
    size_t N = 0, N2 = 0, M = 0;

    double nCont = args->nCont;
    size_t nStruct = args->nStruct;

    pargs_show(stdout, args);
    printf("Press <Enter> to continue\n");
    getchar();

    uint8_t * L = readl(lFile, &M);
    uint8_t lmin = 255; uint8_t lmax = 0;
    for(size_t kk = 0; kk<M; kk++)
    {
        if(lmin > L[kk])
            lmin = L[kk];
        if(lmax < L[kk])
            lmax = L[kk];
    }
    printf("    L has %zu elements, min: %u, max %u\n", M, lmin, lmax);

    double * H = readh(hFile, &N2);
    N = sqrt(N2);
    assert(N*N == N2);
    assert(N==M);
    printf("    H has %zu elements, [%zu x %zu]\n", N2, N, N);


    size_t nmax = 0;
    double cmax = -1e99;
    double cmean = mostContacts(H, NULL, N, nStruct, &cmax, &nmax);
    printf(" >> Before any processing:\n");
    printf("    Most contacts for bead %zu : %f\n", nmax, cmax);
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
        N2 = N*N;
    }

    printf(" >> Removing beads where N(i,i) < max(N(i,:))\n");
    clearWeak(H, N);


    if(args->mode == MODE_DEFAULT)
    {
        printf("    Setting diag(H) = 0\n");
        clearDiagonal(H, N, 0);

        printf(" >> Finding scaling factor\n");
        scale(H, L, N, nStruct, nCont);
    }

    if(args->mode == MODE_EQ)
    {
        printf("    EQUAL-MODE preparations\n");
        clearDiagonal(H, N, 1);

        printf(" >> Balancing H\n");
        double berror = balance(H, N);
        printf("    Largest error: %f\n", berror);
        //  writeA(H, N2, "H_balanced.double");

        printf(" >> Finding scaling factor\n");
        scale(H, L, N, nStruct, nCont);
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
