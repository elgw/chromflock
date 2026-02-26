//

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include "cf_util.h"
#include "gzl.h"

#include "txt/con2mflock_progdesc.txt.h"
#include "txt/con2mflock_changelog.txt.h"
#include "txt/con2mflock_examples.txt.h"

#define CON2MFLOCK_VERSION_MAJOR 1
#define CON2MFLOCK_VERSION_MINOR 0
#define CON2MFLOCK_VERSION_PATCH 1

typedef int64_t i64;
typedef uint32_t u32;
typedef uint8_t u8;

typedef struct options {
    i64 resolution;
    char * infile; /* .pairs.gz or .con */
    int pairs_format; // set to 1 if a pairs file

    /* These are 0-indexed, i.e. chr1 is at index 0 */
    i64 * chr_sizes;
    i64 * chr_reads;

    i64 nchr; // number of non-zero entries in chr_sizes and chr_reads

    /* Output files */
    char * label_file;
    char * contact_file;
    u32 * contacts;
    i64 nlines;
    int gen_matrix;
    char * matrix_file;
    int verbose;
    // TODO:
    char * genome; // Either a genome file (json?) or the name of a genome
} opts;

opts * opts_new(void)
{
    opts * s = calloc(1, sizeof(opts));
    s->resolution = 1000000;
    s->chr_sizes = calloc(50, sizeof(i64));
    s->chr_reads = calloc(50, sizeof(i64));
    s->nchr = 23;
    s->label_file = strdup("c2m_labels.npy");
    s->contact_file = strdup("c2m_contacts.npy");
    s->matrix_file = strdup("c2m_matrix.npy");
    s->genome = strdup("T2T");
    s->verbose = 1;
    return s;
}

void opts_free(opts * s)
{
    free(s->chr_sizes);
    free(s->chr_reads);
    free(s->infile);
    free(s->label_file);
    free(s->contact_file);
    free(s->contacts);
    free(s->matrix_file);
    free(s->genome);
    free(s);

    return;
}

i64 parse_chr_id(const char * S)
{
    const char * P = S;
    while(*P != 'r' && *P != '\0'){ P++; }
    P++;

    if(P[0] == 'X')
    {
        return 22;
    }
    if(P[0] == 'Y')
    {
        return 23;
    }
    return atol(P) - 1;
}


typedef struct {
    char * name; // "--help"
    int has_arg; // no_argument=0, required_argument=1, option_argument=2
    int * flag; //
    int val; // 'h'
    char * help; // help message
} cmdopt;

static void show_help(char * progname, cmdopt * options)
{
    printf("%s\n", __con2mflock_progdesc_txt);
    printf("Usage: %s [options]\n", progname);
    printf("\n");
    printf("These are the options:\n");
    int i = 0;
    while(options[i].name != NULL)
    {
        printf("  --%s", options[i].name);
        if(options[i].has_arg)
        {
            printf(" arg");
        }
        printf("\n");
        printf("\t%s\n", options[i].help);
        printf("\n");
        i++;
    }
}

void parse_command_line(int argc, char ** argv, opts * s)
{
    cmdopt options[] = {
        {"resolution", required_argument,  NULL, 'r',
         "The binning resolution, i.e. number of basepairs per bin"},
        {"con",        required_argument,  NULL, 'c',
         ".con file to parse"},
        {"pairs",      required_argument,  NULL, 'p',
         ".pairs(.gz) file to parse"},
        {"help",       no_argument,        NULL, 'h',
         "Show this help message"},
        {"matrix",     no_argument,        NULL, 'm',
         "Also generate a contact matrix and write to disk"},
        {"mfile",      required_argument,  NULL, 'M',
         "Set the file name for the matrix output file, also enables --matrix"},
        {"genome",     required_argument,  NULL, 'g',
         "Default: 'T2T', other possible options: 'hg19'. For other genome please add them to the source."},
        {"verbose",    required_argument,  NULL, 'v',
         "Set verbose level"},
        {"version",    no_argument,        NULL, 'V',
         "Show version history"},
        {NULL, 0, NULL, 0, NULL}
    };

    int nopt = 0;
    while(options[nopt].name != NULL)
    {
        nopt++;
    }

    struct option * longopts = calloc(nopt+1, sizeof(struct option));
    for(i64 kk = 0; kk < nopt; kk++)
    {
        longopts[kk].name = options[kk].name;
        longopts[kk].has_arg = options[kk].has_arg;
        longopts[kk].flag = options[kk].flag;
        longopts[kk].val = options[kk].val;
    }

    int ch;
    while((ch = getopt_long(argc, argv,
                            "r:c:hmMp:g:v:V",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'g':
            free(s->genome);
            s->genome = strdup(optarg);
            break;
        case 'r':
            s->resolution = atol(optarg);
            break;
        case 'c':
            free(s->infile);
            s->infile = strdup(optarg);
            s->pairs_format = 0;
            break;
        case 'h':
            show_help(argv[0], options);
            exit(EXIT_SUCCESS);
            break;
        case 'm':
            s->gen_matrix = 1;
            break;
        case 'M':
            s->gen_matrix = 1;
            free(s->matrix_file);
            s->matrix_file = strdup(optarg);
            break;
        case 'p':
            free(s->infile);
            s->infile = strdup(optarg);
            s->pairs_format = 1;
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case 'V':
            printf("%s", __con2mflock_changelog_txt);
            exit(EXIT_SUCCESS);
        }
    }

    if(s->infile == NULL)
    {
        fprintf(stderr, "Error: Please provide --pairs or --confile (none given)\n\n");
        show_help(argv[0], options);
        exit(EXIT_FAILURE);
    }
    free(longopts);
    return;
}

static int
parse_con_line(char * L, i64 * chr1, i64 * pos1, i64 * chr2, i64 * pos2)
{
    i64 p = 0;
    i64 ipos = 1;
    i64 pos[6];
    pos[0] = 0;

    while( L[p] != '\n' && L[p] != '\0')
    {
        if(L[p] == ',')
        {
            pos[ipos++] = p+1;
            L[p] = '\0';
        }
        p++;
    }

    if(ipos != 5)
    {
        return 1;
    }
    *chr1 = parse_chr_id(L+pos[0]);
    *pos1 = atol(L+pos[1]);
    *chr2 = parse_chr_id(L+pos[2]);
    *pos2 = atol(L+pos[3]);
    return 0;
}

static int
parse_pairs_line(char * L, i64 * chr1, i64 * pos1, i64 * chr2, i64 * pos2)
{
    //printf("Line = %s\n", L);
    if(L[0] == '#')
    {
        return 1;
    }


    i64 p = 0;
    i64 ipos = 1;
    i64 pos[8];
    pos[0] = 0;

    while( L[p] != '\n' && L[p] != '\0')
    {
        if(L[p] == '\t')
        {
            pos[ipos++] = p+1;
            L[p] = '\0';
        }
        p++;
    }

    *chr1 = parse_chr_id(L+pos[1]);
    *pos1 = atol(L+pos[2]);
    *chr2 = parse_chr_id(L+pos[3]);
    *pos2 = atol(L+pos[4]);

    //printf("%ld, %ld, %ld, %ld\n", chr1[0], pos1[0], chr2[0], pos2[0]);
    return 0;
}

// GRCh38.p14
// https://www.ncbi.nlm.nih.gov/grc/human/data
i64 chr_size_hg38[] =
    {248956422, 242193529, 198295559, 190214555,
     181538259, 170805979, 159345973, 145138636,
     138394717, 133797422, 135086622, 133275309,
     114364328, 107043718, 101991189,  90338345,
     83257441,   80373285,  58617616,  64444167,
     46709983,   50818468, 156040895,  57227415}; // 21, 22, X, Y

// T2T-CHM13v1.1
// https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.3/
i64 chr_size_CHM13[] =
    {248387328, 242696752, 201105948, 193574945, // 1-4
     182045439, 172126628, 160567428, 146259331, // 5-8
     150617247, 134758134, 135127769, 133324548, // 9-12
     113566686, 101161492,  99753195,  96330374, // 13-16
     84276897,  80542538,  61707364,  66210255, // 17-20
     45090682,  51324926, 154259566};           // 21, 22, X

static void get_chr_size(opts * s)
{
    assert(s->genome != NULL);
    int chr_size_fixed = 0;

    if(strcmp(s->genome, "hg19") == 0)
    {
        for(i64 kk = 0 ; kk < 23; kk++)
        {
            s->chr_sizes[kk] = chr_size_hg38[kk];
        }
        chr_size_fixed = 1;
    }
    if(strcmp(s->genome, "T2T") == 0)
    {
        for(i64 kk = 0 ; kk < 23; kk++)
        {
            s->chr_sizes[kk] = chr_size_CHM13[kk];
        }
        chr_size_fixed = 1;
    }

    if(chr_size_fixed == 0)
    {
        printf("There is no chromosome sizes available for %s\n", s->genome);
        exit(EXIT_FAILURE);
    }
    return;
}

static i64 get_dataset_size(const opts * s)
{
    // Side effects:
    // Sets s->n_lines

    int chr_size_fixed = 1;
    if(s->verbose > 0)
    {
        printf("Checking number of reads\n");
    }


    gzl_state * gzl = gzl_open(s->infile, 1024);

    if(gzl == NULL)
    {
        printf("Failed to read %s\n", s->infile);
        exit(EXIT_FAILURE);
    }

    char * L = NULL;
    i64 nlines = 0;
    int gzl_error = 0;


    while( (L = (char *) gzl_get_line(gzl, &gzl_error) ) )
    {
        nlines++;
        i64 chr1, chr2, pos1, pos2;


        if(s->pairs_format)
        {
            if(L[0] == '#')
            {
                continue;
            }
            parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2);
        } else {
            if(parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2))
            {
                printf("Unable to parse a line #%ld\n", nlines+1);
                continue;
            }
        }
        if( (chr1 >= s->nchr) || (chr2 >= s->nchr))
        {
            printf("Parsed something weird: chr1: %ld, chr2: %ld\n", chr1, chr2);
            printf("Line %ld: '%s'\n", nlines, L);
            exit(EXIT_FAILURE);
        }
        s->chr_reads[chr1]++;
        s->chr_reads[chr2]++;
        if(chr_size_fixed)
        {
            if(s->chr_sizes[chr1] < pos1)
            {
                fprintf(stderr, "ERROR: Got (chrid=%ld, pos=%ld) but that chr is only %ld large\n",
                        chr1, pos1, s->chr_sizes[chr1]);
                exit(EXIT_FAILURE);
            }
            if(s->chr_sizes[chr2] < pos2)
            {
                fprintf(stderr, "ERROR: Got (chrid=%ld, pos=%ld) but that chr is only %ld large\n",
                        chr2, pos2, s->chr_sizes[chr2]);
                exit(EXIT_FAILURE);
            }
        } else {
            #if 0
            if(s->chr_sizes[chr1] < pos1)
            {
                s->chr_sizes[chr1] = pos1;
            }
            if(s->chr_sizes[chr2] < pos2)
            {
                s->chr_sizes[chr2] = pos2;
            }

            if(chr1 > s->nchr)
            {
                s->nchr = chr1+1;
            }
            if(chr2 > s->nchr)
            {
                s->nchr = chr2+1;
            }
            #endif
        }
        nlines++;
    }

    if(s->nchr == 0)
    {
        fprintf(stderr, "\nERROR: Could not get any data from the file\n");
        exit(EXIT_FAILURE);
    }


    printf("Read %ld lines\n", nlines);
    printf(" #,      Size,    Reads,    Bins\n");

    for(i64 kk = 0; kk < s->nchr; kk++)
    {
        printf("%2ld, %9ld, %8ld, %7ld\n", kk+1,
               s->chr_sizes[kk], s->chr_reads[kk],
               (s->resolution + s->chr_sizes[kk])/s->resolution);
    }
    gzl_destroy(gzl);
    free(L);
    return nlines;
}

/* Write the label array as uint8
 * chromflock is happy if we leave 0 unused
 */

static void
write_labels(const opts * s)
{
    printf("Writing labels to %s\n", s->label_file);
    // number of labels
    size_t nbin = 0;
    for(i64 cc = 0; cc < s->nchr; cc++)
    {
        i64 n = (s->chr_sizes[cc] + s->resolution)/s->resolution;
        nbin += n;
    }

    // put in array
    u8 * labels = malloc(nbin*sizeof(u8));
    size_t idx = 0;
    for(i64 cc = 0; cc < s->nchr; cc++)
    {
        u8 chr = cc+1;
        i64 n = (s->chr_sizes[cc] + s->resolution)/s->resolution;
        for(i64 nn = 0 ; nn < n; nn++)
        {
            labels[idx++] = chr;
        }
    }

    // write to disk
    if(write_bead_labels(s->label_file, labels, nbin))
    {
        printf("Failed to write to %s\n", s->label_file);
        exit(EXIT_FAILURE);
    }

    free(labels);
    return;
}

int cmp_u32_pair(const void * _A, const void * _B)
{
    u32 * A = (u32*) _A;
    u32 * B = (u32*) _B;

    if(A[0] < B[0])
    {
        return -1;
    }

    if(A[0] > B[0])
    {
        return 1;
    }

    if(A[1] < B[1])
    {
        return -1;
    }

    if(A[1] > B[1])
    {
        return 1;
    }

    return 0;

}

static void
write_contacts(opts * s)
{
    assert(s->nchr > 0);
    assert(s->nlines > 0);

    if(s->verbose > 0)
    {
        printf("Writing contacts to %s\n", s->contact_file);
    }

    gzl_state * gzl = gzl_open(s->infile, 1024);

    if(gzl == NULL)
    {
        gzl_destroy(gzl);
        printf("Failed to read %s\n", s->infile);
        exit(EXIT_FAILURE);
    }

    s->contacts = calloc(s->nlines*2, sizeof(u32));

    size_t nchr = s->nchr;
    if(nchr >= 50)
    {
        exit(EXIT_FAILURE);
    }
    i64 * bin_start = calloc(nchr+1, sizeof(i64));

    for(i64 kk = 1; kk < s->nchr; kk++)
    {
        bin_start[kk] = bin_start[kk-1] + (s->chr_sizes[kk-1]+s->resolution)/s->resolution;
    }

    if(s->verbose > 1)
    {
        for(i64 kk = 0; kk < s->nchr; kk++)
        {
            printf("chr #%ld starts at %ld\n", kk+1, bin_start[kk]);
        }
    }

    i64 line = 0;
    char * L = NULL;
    int gzl_error;

    while( (L = (char*) gzl_get_line(gzl, &gzl_error) ) )
    {
        line++;
        i64 chr1, chr2, pos1, pos2;
        if(s->pairs_format)
        {
            if(L[0] == '#')
            {
                continue;
            }
            parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2);

        } else {
            if(parse_con_line(L, &chr1, &pos1, &chr2, &pos2))
            {
                printf("Unable to parse a line #%ld\n", line+1);
                exit(EXIT_FAILURE);
            }
        }
        if((chr1 > s->nchr) || (chr2 > s->nchr))
        {
            printf("Input error on line %ld\n", line);
            printf("chr1=%ld, chr2=%ld\n", chr1, chr2);
            assert(0);
        }


        i64 bead1 = bin_start[chr1] + pos1/s->resolution;
        i64 bead2 = bin_start[chr2] + pos2/s->resolution;

        s->contacts[2*line + 0] = bead1;
        s->contacts[2*line + 1] = bead2;

    }
    gzl_destroy(gzl);

    printf("Sorting contacts\n");
    qsort(s->contacts, s->nlines, 2*sizeof(u32), cmp_u32_pair);

    printf("Removing duplicates\n");
    i64 wpos = 0;
    for(i64 kk = 1; kk < s->nlines; kk++)
    {
        if(cmp_u32_pair(s->contacts + 2*kk, s->contacts + 2*kk-2) != 0)
        {
            memcpy(s->contacts+2*wpos, s->contacts+2*kk, 2*sizeof(u32));
            wpos++;
        }
    }
    printf("Keeping %ld / %ld\n", wpos, s->nlines);

    printf("Writing contacts to %s\n", s->contact_file);

    if(write_u32(s->contact_file, s->contacts, 2, wpos))
    {
        printf("Could not write to %s\n", s->contact_file);
        exit(EXIT_FAILURE);
    }

    free(L);
    free(bin_start);
    return;
}

static void write_matrix(opts * s)
{
    printf("Writing contact map to %s\n", s->matrix_file);

    size_t nchr = s->nchr;
    if(nchr >= 50)
    {
        exit(EXIT_FAILURE);
    }
    i64 * bin_start = calloc(nchr, sizeof(i64));

    for(i64 kk = 1; kk < s->nchr; kk++)
    {
        bin_start[kk] = bin_start[kk-1] + (s->chr_sizes[kk-1]+s->resolution)/s->resolution;
    }

    i64 N = 0;
    for(i64 kk = 0; kk < s->nchr; kk++)
    {
        N+= (s->chr_sizes[kk]+s->resolution)/s->resolution;
    }

    printf("Matrix size: %ld x %ld\n", N, N);

    u32 * M = calloc(N*N, sizeof(u32));

    gzl_state * gzl = gzl_open(s->infile, 1024);

    if(gzl == NULL)
    {
        printf("Failed to read %s\n", s->infile);
        exit(EXIT_FAILURE);
    }

    i64 line = 0;
    char * L = NULL;
    int gzl_error;
    while( (L = (char *) gzl_get_line(gzl, &gzl_error) ) )
    {
        line++;
        i64 chr1, chr2, pos1, pos2;
        if(s->pairs_format)
        {
            if(L[0] == '#')
            {
                continue;
            }
            parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2);
        } else {
            if(parse_con_line(L, &chr1, &pos1, &chr2, &pos2))
            {
                printf("Unable to parse a line #%ld\n", line);
                exit(EXIT_FAILURE);
            }
        }
        i64 bead1 = bin_start[chr1] + pos1/s->resolution;
        i64 bead2 = bin_start[chr2] + pos2/s->resolution;
        M[bead1 + N*bead2] += 1;

    }
    gzl_destroy(gzl);

    // Make symmetric
    // Note: We don't touch the diagonal

    for(i64 kk = 0; kk < N; kk++)
    {
        for(i64 ll = kk +1; ll < N; ll++)
        {
            i64 reads = M[kk + N*ll] + M[ll + N*kk];
            M[kk + N*ll] = reads;
            M[ll + N*kk] = reads;
        }
    }

    printf("Writing to %s\n", s->matrix_file);
    if(write_u32(s->matrix_file, M, N, N))
    {
        printf("Failed to write to %s!\n", s->matrix_file);
        exit(EXIT_FAILURE);
    }

    free(M);
    free(L);
    free(bin_start);
    return;
}


int con2mflock(int argc, char ** argv)
{
    opts * s = opts_new();

    parse_command_line(argc, argv, s);

    get_chr_size(s);

    // Check how many lines of input
    s->nlines = get_dataset_size(s);

    /* with fixed chr sizes we have not even touched the input
       data so while this works write_contacts might fail leaving us
       with an incomplete set of output files. */
    write_labels(s);

    /* - parses the input data for contacts
     * - sorts the contacts
     * - removes duplicates (the lower the resolution, the more duplicates)
     * - writes in a format that chromflock can read
     */
    write_contacts(s);

    /* Generate a dense contact matrix, not advised for high resolutions */
    if(s->gen_matrix)
    {
        write_matrix(s);
    }

    opts_free(s);
    printf("done\n");
    return EXIT_SUCCESS;
}
