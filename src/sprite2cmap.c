#include "sprite2cmap.h"
#include "balance.h"

#define UNUSED(...) (void)(__VA_ARGS__)

typedef struct {
    int balance; // Matrix Balancing, yes/no
    size_t binSize; // Size of each bin
    size_t maxClusterSize;
    size_t minClusterSize;
    int printBins;
    int ok; // skip all questions
    int fileType;
    int chr; // if > 0 only extract that chr
    int verbose;
    char * bFile;
    char * cFile; // .cluster file (input)
    char * cFileBase; // basename of cFile
    char * logFile;
    char * oFile; // output file, heatmat
    char * coordFile; // output file, coordinates
    char * lFile; // uint8 file with chromosome labels per bead
} pset;


static void usage(char *);
static void pset_show(pset *, FILE *);
static pset * pset_alloc(void);
static void pset_free(pset * p);


// hg19 chromosomes sizes from
// https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
size_t chrSize[] = {
    0, // DUMMY
    249250621, // chr1
    243199373, // chr2
    198022430, // chr3
    191154276, // chr4
    180915260, // chr5
    171115067, // chr6
    159138663, // chr7
    146364022, // chr8
    141213431, // chr9
    135534747, // chr10
    135006516, // chr11
    133851895, // chr12
    115169878, // chr13
    107349540, // chr14
    102531392, // chr15
    90354753 , // chr16
    81195210 , // chr17
    78077248 , // chr18
    59128983 , // chr19
    63025520 , // chr20
    48129895 , // chr21
    51304566 , // chr22
    155270560, // chrX  (23)
    59373566 , // chrY  (24)
};


// This is not really "error" handling just a direction of the flow before it crashes the program
// would be easier to read if it was just a function call ...
jmp_buf ebuf;

int eline = 0;
size_t emem = 0;

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

#if 0
#define malloc(x) tmalloc(x, __FILE__, __LINE__)
#define calloc(x, y) tcalloc(x, y, __FILE__, __LINE__)
static void * tcalloc(size_t nmemb, size_t size, char * file, int line)
{
    UNUSED(file);
    void * A = (calloc)(nmemb, size); // The () around calloc says that the macro shouldn't be used.
    if(A == 0)
    {
        emem = nmemb*size;
        eline = line;
        longjmp(ebuf, 1);
    }
    return A;
}

static void * tmalloc(size_t nbytes, char * file, int line)
{
    UNUSED(file);
    void * A = (malloc)(nbytes);
    if(A == 0)
    {
        emem = nbytes;
        eline = line;
        longjmp(ebuf, 1);
    }
    return A;
}
#endif

static void writeM(double * M, // data
                   size_t N, // number of elements
                   char * outName)
{
    FILE * fout = fopen(outName, "w");
    fwrite(M, sizeof(double), N, fout);
    fclose(fout);
}

static int parsechr(char * str, int * chr)
{
    if( (strlen(str) < 3) || (strlen(str) > 5))
    {
        return 0;
    }

    if(str[3] == 'X')
    {
        chr[0] = 23;
        return(1);
    }

    if(str[3] == 'Y')
    {
        chr[0] = 24;
        return(1);
    }

    if(str[3] < '0' || str[3] > '9')
    {
        return(0);
    }

    chr[0] = atoi(str+3);
    return(1);
}

static int parsepos(char * str, size_t * pos)
{
    pos[0] = atol(str);
    return 1;
}

static int parseSpriteWord(char * str, int * chr, size_t * pos)
{
    //  printf("to parse: %s\n", str);
    size_t len = strlen(str);
    char * posStr = NULL;
    int ok = 0;

    for(size_t kk = 0 ; kk<len; kk++)
    {
        if(str[kk] == ':')
        {
            str[kk] = '\0';
            posStr = str + kk + 1;
            ok = 1;
        }
    }

    if(ok == 0)
    {
        printf("Warning: can't parse string '%s'\n", str);
        return(0);
    }

    ok = parsechr(str, chr);
    if(ok == 0)
    {
        return 0;
    }

    ok = parsepos(posStr, pos);
    if(ok == 0)
    {
        return 0;
    }
    return 1;
}

static size_t getBin(int * S, int * E, size_t binSize, int chr, size_t pos)
{
    UNUSED(E);
    size_t bin = S[chr] + floor(pos/binSize);
    assert(bin <= (size_t) E[chr]);

    return bin;
}

static int size_t_cmp(const void * a, const void * b)
{
    return ((size_t *) a)[0] > ((size_t *) b)[0];
}

static void printBins(FILE * fd, int * S, int *E,
                      size_t binSize,
                      size_t nMemb, int * chr, size_t * pos)
{
    /* Print the bins that the elements of the clusters would go to.
     * Only do this if the bins are unique */

    size_t * bpos = malloc(nMemb * sizeof(size_t));
    assert(bpos != 0);

    for(size_t kk = 0; kk<nMemb; kk++)
    {
        bpos[kk] = getBin(S, E, binSize, chr[kk], pos[kk]);
    }
    qsort(bpos, nMemb, sizeof(size_t), size_t_cmp);

    int unique = 1;
    for(size_t kk = 1; kk<nMemb; kk++)
    {
        if(bpos[kk] == bpos[kk-1])
        {
            unique = 0;
        }
    }
    if(unique)
    {
        fprintf(fd, "%zu", bpos[0]);
        for(size_t kk = 1; kk<nMemb; kk++)
        {
            fprintf(fd, ", %zu", bpos[kk]);
        }

        fprintf(fd, "\n");
    }
    free(bpos);

}

static void addToMap(double * M, size_t N,  // M: map of size NxN
                     int * S, int * E, // Start and end positions of the chromosomes
                     size_t binSize,
                     size_t nMemb, int * chr, size_t * pos) // Size of cluster and positions
{

    /* See the section 'Generating pairwise contacts and heatmaps from sprite data'
     * in the paper */

    const double weight = 2.0/(double) nMemb;

    for(size_t kk = 0; kk<nMemb; kk++)
    {
        if(chrSize[chr[kk]] > 0)
        {
            size_t pos1 = getBin(S, E, binSize, chr[kk], pos[kk]);
            for(size_t  ll = kk+1; ll<nMemb; ll++)
            {
                if(chrSize[chr[ll]] > 0)
                {
                    size_t pos2 = getBin(S, E, binSize, chr[ll], pos[ll]);
                    assert(pos1<N);
                    assert(pos2<N);

                    M[pos1 + N*pos2]+=weight;
                    M[pos2 + N*pos1]+=weight;
                }
            }
        }
    }
    return;
}


static int argparsing(pset * p, int argc, char ** argv)
{
    struct option longopts[] = {
        {"version",     no_argument,       NULL,   'i' },
        { "help",         no_argument,       NULL,   'h' },
        {"defaults",      no_argument,       NULL,  'd' },
        {"balance",      no_argument,       NULL,  'b' },
        { "ok",           no_argument, NULL, 'k'},
        // Data
        { "binSize", required_argument, NULL, 'z'},
        { "maxSize",        required_argument, NULL,   'c' },
        { "minSize",            required_argument, NULL, 'C' },
        { "sFile",        required_argument, NULL,   's' },
        { "oFile",        required_argument, NULL,   'o' },
        { "chr",          required_argument, NULL, 'S' },
        { "printBins",    no_argument, NULL, 'B'},
        { "hic",          no_argument, NULL, 'H'},
        { "verbose",      required_argument, NULL, 'v'},
        { NULL, 0, NULL, 0}
    };

    int ch;
    while((ch = getopt_long(argc, argv, "z:bBcdhHic:C:o:s:S:v:", longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'i':
            printf("Build date: %s, %s\n", __DATE__, __TIME__);
#ifdef GIT_VERSION
            printf("GIT HASH: %s\n", GIT_VERSION);
#endif
#ifdef CC_VERSION
            printf("Compiler: %s\n", CC_VERSION);
#endif
            exit(0);
        case 'd':
            printf("Defaults:\n");
            pset_show(p, stdout);
            exit(0);
            break;
        case 'v':
            p->verbose = atoi(optarg);
            break;
        case 'B':
            p->printBins = 1;
            break;
        case 'S':
            p->chr = atoi(optarg);
            break;
        case 'c':
            p->maxClusterSize = atol(optarg);
            break;
        case 'C':
            p->minClusterSize = atol(optarg);
            break;
        case 's':
            if(p->cFile != NULL){free(p->cFile);}
            p->cFile = malloc(strlen(optarg)+1);
            assert(p->cFile != NULL);
            strcpy(p->cFile, optarg);
            char * temp = strdup(p->cFile);
            char * temp2 = basename(temp); // don't free this
            p->cFileBase = strdup(temp2);
            free(temp);
            break;
        case 'o':
            if(p->oFile != NULL){free(p->oFile);}
            p->oFile = malloc(strlen(optarg)+1);
            assert(p->oFile != NULL);
            strcpy(p->oFile, optarg);
            break;
        case 'b':
            p->balance = 1;
            break;
        case 'k':
            p->ok = 1;
            break;
        case 'H':
            p->fileType = 1;
            break;
        case 'z':
            p->binSize = atol(optarg);
            break;
        default:
            printf("Unknown command line option\n");
            return(1);
        }
    }

    if(p->cFile == NULL)
    {
        printf("--sFile not specified\n");
        return -1;
    }

    if(p->printBins == 1)
    {
        p->bFile = malloc((strlen(p->cFileBase)+100)*sizeof(char));
        assert(p->bFile != NULL);
        sprintf(p->bFile, "%s.bins", p->cFileBase);
    }

    if(p->oFile == NULL)
    {
        p->oFile = malloc((strlen(p->cFileBase)+100)*sizeof(char));
        assert(p->oFile != NULL);
        sprintf(p->oFile, "%s.heatmap.double", p->cFileBase);

        p->coordFile = malloc((strlen(p->cFileBase)+100)*sizeof(char));
        assert(p->coordFile != NULL);
        sprintf(p->coordFile, "%s.heatmap.coords.tsv", p->cFileBase);
    }

    if(p->lFile == NULL)
    {
        p->lFile = malloc((strlen(p->cFileBase)+100)*sizeof(char));
        assert(p->lFile != NULL);
        sprintf(p->lFile, "%s.heatmap.L.uint8", p->cFileBase);
    }

    if(p->logFile == NULL)
    {
        p->logFile = malloc(100+strlen(p->oFile));
        sprintf(p->logFile, "%s.log.txt", p->oFile);
    }

    if(p->chr > 0)
    {
        for(int kk = 0; kk<25; kk++)
        {
            if(kk != p->chr)
            {
                chrSize[kk] = 0;
            }
        }
    }

    return 0;
}

static void usage(char * cmdName)
{
    printf("%s is used to convert SPRITE (.cluster) files to heatmaps. It can also convert some other files to heat maps.\n", cmdName);
    printf("Usage:\n");
    printf("%s ...\n", cmdName);
    printf("Required arguments:\n");
    printf(" --sFile <file>\n\tsprite file (.cluster)\n");
    printf("Optional arguments\n");
    printf(" --oFile <file>\n\toutput file\n");
    printf(" --binSize N\n\tsize of each bin\n");
    printf(" --balance\n\tenable matrix balancing\n");
    printf(" --maxSize s\n\tmax cluster size\n");
    printf(" --minSize s\n\tmin cluster size\n");
    printf(" --printBins\n\tGenerate a .bins file where genomic coordinates are converted to bin number.\n");
    printf(" --hic\n\tParse FAHiC-files\n");
    printf(" --chr c\n\tOnly extract chr c\n");
    printf(" --verbose N\n\tSet verbose level to N, default=0\n");
    printf(" --version\n\tShow version\n");
    return;
}

static pset * pset_alloc(void)
{
    pset * p = calloc(1, sizeof(pset));
    p->balance = 0;
    p->binSize = 1000000;
    p->maxClusterSize = 1000;
    p->minClusterSize = 2;
    p->ok = 0;
    p->cFile = NULL;
    p->printBins = 0;
    p->bFile = NULL;
    p->cFileBase = NULL;
    p->logFile = NULL;
    p->oFile = NULL;
    p->coordFile = NULL;
    p->lFile = NULL;
    p->fileType = 0;
    p->chr = 0;
    return p;
}

static void pset_show(pset * p, FILE * fout)
{
    fprintf(fout, "Bin size: %zu\n", p->binSize);
    fprintf(fout, "Balancing: %d\n", p->balance);
    fprintf(fout, "Min cluster Size: %zu\n", p->minClusterSize);
    fprintf(fout, "Max cluster Size: %zu\n", p->maxClusterSize);
    fprintf(fout, "Input file: %s\n", p->cFile);
    fprintf(fout, "Output file: %s\n", p->oFile);
    fprintf(fout, "Log file: %s\n", p->logFile);
    fprintf(fout, "Label file: %s\n", p->lFile);
    fprintf(fout, "Coordinate file: %s\n", p->coordFile);
    fprintf(fout, "File type: %d\n", p->fileType);
}

static void pset_free(pset * p)
{
    if(p->bFile != NULL)
    {
        free(p->bFile);
        p->bFile = NULL;
    }
    if(p->cFile != NULL)
    {
        free(p->cFile);
        p->cFile = NULL;
    }
    if(p->cFileBase != NULL)
    {
        free(p->cFileBase);
        p->cFileBase = NULL;
    }
    if(p->oFile != NULL)
    {
        free(p->oFile);
        p->oFile = NULL;
    }
    if(p->coordFile != NULL)
    {
        free(p->coordFile);
        p->coordFile = NULL;
    }
    if(p->lFile != NULL)
    {
        free(p->lFile);
        p->lFile = NULL;
    }

    if(p->logFile != NULL)
    {
        free(p->logFile);
        p->logFile = NULL;
    }


    free(p);
    return;
}

static void writeCoords(const int * S, const int * E, size_t N, const pset * p)
{

    printf("Writing coordinates to %s\n", p->coordFile);
    printf(" and chromosome labels to %s\n", p->lFile);

    FILE * cOut = fopen(p->coordFile, "w");
    assert(cOut != NULL);
    FILE * lOut = fopen(p->lFile, "w");
    assert(lOut != NULL);

    fprintf(cOut, "bin\tchr\tstart\tend\t\n");
    for(size_t kk = 0; kk < N; kk++)
    {
        size_t chr = 0;
        for(size_t ll = 1; ll<25; ll++)
        {
            if( ((int) kk >= S[ll]) && ((int) kk <= E[ll]))
            {
                chr = ll;
            }
        }

        assert(chr>0); assert(chr<25);

        size_t from = (kk-S[chr])*p->binSize;
        size_t to = from+p->binSize-1;
        fprintf(cOut, "%zu\t%zu\t%zu\t%zu\n", kk, chr, from, to);
        uint8_t chr_u8 = chr;

        fwrite(&chr_u8, sizeof(uint8_t), 1, lOut);
    }

    fclose(cOut);
    fclose(lOut);
}

static void writeSettings(pset * p, int argc, char ** argv)
{
    printf("Writing settings to %s\n", p->logFile);
    FILE * f = fopen(p->logFile, "w");
    assert(f != NULL);


    fprintf(f, "# Build date: %s, %s\n", __DATE__, __TIME__);
#ifdef GIT_VERSION
    fprintf(f, "# GIT HASH: %s\n", GIT_VERSION);
#endif
#ifdef CC_VERSION
    fprintf(f, "# Compiler: %s\n", CC_VERSION);
#endif
    fprintf(f, "# Command line:\n");
    for(int kk = 0 ; kk<argc; kk++)
    {
        fprintf(f, "%s ", argv[kk]);
    }
    fprintf(f, "\n");
    fclose(f);
}

static int parseFAHIC(char * buf, // Line to parse
                      char * buf2, // Buffer
                      size_t maxClusters,
                      // Return:
                      int * chr, // Chromosomes
                      size_t * pos, // Position
                      size_t * N) // Number of hits
{
    UNUSED(buf2);
    UNUSED(maxClusters); // All clusters are expected to have two members

    if(buf[0] == '#')
    {
        N[0] = 0;
        return 1;
    }

    char * found = buf;
    char * next = buf;

    found = strsep(&next, "\t");

    found = strsep(&next, "\t");
    parsechr(found, chr);
    found = strsep(&next, "\t");
    parsepos(found, pos);

    found = strsep(&next, "\t");
    parsechr(found, chr+1);
    found = strsep(&next, "\t");
    parsepos(found, pos+1);

    //    printf("%d/%zu %d/%zu\n", chr[0], pos[0], chr[1], pos[1]);

    N[0] = 2;

    return 1;
}

static int parseSprite(char * buf, // Line to parse
                       char * buf2, // Buffer
                       size_t maxClusters,
                       // Return:
                       int * chr, // Chromosomes
                       size_t * pos, // Position
                       size_t * N) // Number of hits
{
    size_t nMemb = 0;
    int firstWord = 1;
    char * found = buf;
    char * next = buf;

    while( (found = strsep(&next,"\t")) != NULL && (nMemb <= maxClusters))
    {
        if(firstWord==0)
        {
            strncpy(buf2, found, 100);
            if(parseSpriteWord(buf2, chr+nMemb, pos+nMemb) == 1)
            {
                nMemb++;
            }
        }
        firstWord = 0;
    }

    N[0] = nMemb;
    return 1;
}

static int setStartEnd(pset * p, int * S, int * E)
{
    /* Sets the start and end bin of each chromosome into the arrays S and E using
     * chrSize and p->binSize
     * Returns the total number of bins.
     */

    size_t N = 0;

    S[0] = 0; E[0] = 0;

    for(int kk = 1; kk<25; kk++)
    {
        if(chrSize[kk-1] > 0)
        {
            S[kk] = E[kk-1] + 1;
        } else {
            S[kk] = E[kk-1];
        }

        E[kk] = S[kk] + floor((double) chrSize[kk]/(double) p->binSize);

        if(chrSize[kk] > 0)
        {
            N += E[kk]-S[kk] + 1;
        }

        if(p->verbose > 1)
        {
            if(chrSize[kk] > 0)
            {
                printf("chr#%d, S=%d, E=%d\n", kk, S[kk], E[kk]);
            }
        }
    }

    return N;
}


static int convert(int argc, char ** argv)
{

    pset *p = pset_alloc();

    if(argparsing(p, argc, argv))
    {
        usage(argv[0]);
        pset_free(p);
        return 1;
    }

    printf("Parameters:\n");
    pset_show(p, stdout);
    if(p->ok == 0)
    {
        printf("Press <enter> to continue\n");
        getchar();
    }

    FILE * fin = fopen(p->cFile, "r");
    if(fin == NULL)
    {
        printf("Unable to open %s\n", p->cFile);
        return(-1);
    }

    writeSettings(p, argc, argv);

    /* Set up bins in matrix M as well as a mapping from
     * (chr, pos) -> binNumber
     */
    int * S = malloc(25*sizeof(int));
    int * E = malloc(25*sizeof(int));
    size_t N = setStartEnd(p, S, E);
    // N: total number of bins

    writeCoords(S, E, N, p);

    printf("Creating a %zux%zu contact map\n", N, N);
    double * M = calloc(N*N, sizeof(double));

    /*
     * To get largest cluster size:
     awk 'length > max_length { max_length = length; longest_line = $0 } END { print longest_line }' file | wc
    */
    size_t bufMaxSize = 2000*2048;
    char * buf = malloc(bufMaxSize*sizeof(char));
    char * buf2 = malloc(bufMaxSize*sizeof(char));
    size_t maxClusters = bufMaxSize/20; // Max number of members per cluster

    assert(bufMaxSize > p->maxClusterSize*20);

    // Populate these arrays per read line:
    int * chr = malloc(maxClusters*sizeof(int));
    size_t * pos = malloc(maxClusters*sizeof(size_t));

    ssize_t nread;

    printf("Parsing %s\n", p->cFile);
    size_t nTooBig = 0;
    size_t nTooSmall = 0;
    size_t nLines = 0; // counter for number of lines read
    size_t largestCluster = 0;
    const size_t maxClusterSize = p->maxClusterSize;
    const size_t minClusterSize = p->minClusterSize;

    FILE * bFile = NULL;
    if(p->printBins == 1)
    {
        bFile = fopen(p->bFile, "w");
    }

    while((nread = getline(&buf, &bufMaxSize, fin)) != -1)
    {
        nLines++;

        if(nLines % 100000 == 0)
        {
            printf("\rLine: %zu ...", nLines);
            fflush(stdout);
        }

        size_t nMemb = 0;

        int ok = 0;
        if(p->fileType == 0)
        {
            ok = parseSprite(buf, buf2, maxClusters, chr, pos, &nMemb);
        }

        if(p->fileType ==1)
        {
            ok = parseFAHIC(buf, buf2, maxClusters, chr, pos, &nMemb);
        }

        if(ok == 0)
        {
            printf("Can't parse line %zu\n", nLines);
        }

        if(nMemb <= maxClusterSize && nMemb >= minClusterSize)
        {
            if(0)
            {
                for(size_t kk = 0; kk<nMemb; kk++)
                { printf("%d(%zu) ", chr[kk], pos[kk]); }
                printf("\n");
            }

            addToMap(M, N,
                     S, E,
                     p->binSize,
                     nMemb, chr, pos);
            if(p->printBins == 1)
            {
                printBins(bFile, S, E, p->binSize, nMemb, chr, pos);
            }

        } else {
            if(nMemb < minClusterSize)
            {
                nTooSmall++;
            }
            if(nMemb > maxClusterSize)
            {
                nTooBig++;
            }
        }
        if(nMemb > largestCluster)
        {
            largestCluster = nMemb;
        }
    }

    if(bFile != NULL)
    {
        fclose(bFile);
    }

    printf("\rRead %zu lines\n", nLines);
    printf("Excluded clusters: %zu too big, %zu too small\n", nTooBig, nTooSmall);
    printf("Largest cluster: %zu\n", largestCluster);
    if(p->balance == 1)
    {
        printf("Balancing\n");
        double berror = balance(M, N);
        printf("Largest error: %f\n", berror);
    }

    free(chr);
    free(pos);
    fclose(fin);
    free(buf2);
    free(buf);
    free(S);
    free(E);

    writeM(M, N*N, p->oFile);
    free(M);

    printf("Done\n");

    pset_free(p);

    return(0);
}


int sprite2cmap(int argc, char ** argv)
{

    int ret = 0;

    switch(setjmp(ebuf)) // TRY
    {
    case 0:
        ret = convert(argc, argv);
        break;
    case 1: // CATCH
        printf("Failed to allocate %zu MB of memory on line %d\n", emem/100000, eline);
        ret = -1;
        break;
    default:
        printf("Unknown exception thrown on line %d\n", eline);
        ret = -1;
    } // ENDTRY


    return ret;
}
