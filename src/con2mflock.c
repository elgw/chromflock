//

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include "ddict.h"
#include "cf_util.h"
#include "gzl.h"

#include "txt/con2mflock_progdesc.txt.h"
#include "txt/con2mflock_changelog.txt.h"
#include "txt/con2mflock_examples.txt.h"

#define CON2MFLOCK_VERSION_MAJOR 1
#define CON2MFLOCK_VERSION_MINOR 0
#define CON2MFLOCK_VERSION_PATCH 1

typedef struct options {
    i64 resolution;
    int write_backbone;
    char * infile; // .pairs.gz or .con
    int pairs_format; // set to 1 if a pairs file

    // Following arrays are 0-indexed, i.e. chr1 is at index 0
    i64 * chr_size_bp;  // expressed in basepairs
    i64 * chr_size_bin; // as above but expressed in resolution
    i64 * chr_reads; // counts the number of reads that map into each chromosome
    u8 * labels; // chr label per bin

    i64 n_chr; // number of non-zero entries in chr_sizes and chr_reads
    i64 n_bin; // Number of bins/beads in total

    // Output files
    char * label_file;
    char * contact_file;
    char * contact_file_raw;
    u32 * contacts;
    i64 n_reads;
    int gen_matrix;
    char * matrix_file;
    int verbose;
    // A "genome.csv" file to be read by sc_chromflock etc. Provides
    // simply the chromosome names and their size. The line order matters
    // as the first chromsome will have label 1, the 2nd 2 etc.
    char * genome_out_file;

    char * genome; // Either a genome file (json?) or the name of a genome
} opts;

opts * opts_new(void)
{
    opts * s = calloc(1, sizeof(opts));
    s->resolution = 1000000;
    s->chr_size_bp = calloc(50, sizeof(i64));
    s->chr_size_bin = calloc(50, sizeof(i64));
    s->chr_reads = calloc(50, sizeof(i64));
    s->label_file = strdup("c2m_labels.npy");
    s->contact_file = strdup("c2m_contacts.npy");
    s->contact_file_raw = strdup("c2m_contacts_nonbinned.npy");
    s->matrix_file = strdup("c2m_matrix.npy");
    s->genome_out_file = strdup("c2m_genome.csv");
    s->genome = strdup("T2T");
    s->verbose = 1;
    return s;
}

void opts_free(opts * s)
{
    free(s->chr_size_bp);
    free(s->chr_size_bin);
    free(s->labels);
    free(s->chr_reads);
    free(s->infile);
    free(s->label_file);
    free(s->contact_file);
    free(s->contacts);
    free(s->matrix_file);
    free(s->genome);
    free(s->genome_out_file);
    free(s);

    return;
}

// chr1 -> 0, chr2 -> 1 ... chrX->22, chrY -> 23
i64 parse_chr_id(const char * S)
{
    const char * P = S;
    while(*P != 'r' && *P != '\0'){ P++; }
    P++;

    if(P[0] == 'X'){
        return 22;
    }
    if(P[0] == 'Y'){
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

typedef struct {
    u32 chrA;
    u32 posA;
    u32 chrB;
    u32 posB;
} rcontact;

static void
show_help(char * progname, cmdopt * options)
{
    printf("%s\n", con2mflock_progdesc_txt);
    printf("Usage: %s [options]\n", progname);
    printf("\n");
    printf("These are the options:\n");
    int i = 0;
    while(options[i].name != NULL){
        printf("  --%s", options[i].name);
        if(options[i].has_arg){
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
        {"pairs",      required_argument,  NULL, 'p',
         ".pairs(.gz) file to parse. See https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md for the file format specification"},
        {"help",       no_argument,        NULL, 'h',
         "Show this help message"},
        {"matrix",     no_argument,        NULL, 'm',
         "Also generate a contact matrix and write to disk"},
        {"mfile",      required_argument,  NULL, 'M',
         "Set the file name for the matrix output file, also enables --matrix"},
        {"genome",     required_argument,  NULL, 'g',
         "DEPRECIATED. A csv file specifying the genome, i.e. has the columns 'chr_name', and 'chr_size'"},
        {"backbone",   no_argument,        NULL, 'b',
         "Also write 'backbone' contacts, i.e. connect adjacent bins within each chromosome to the contact list (not to the --mfile)"},
        {"verbose",    required_argument,  NULL, 'v',
         "Set verbose level"},
        {"version",    no_argument,        NULL, 'V',
         "Show version history"},
        {NULL, 0, NULL, 0, NULL}
    };

    int nopt = 0;
    while(options[nopt].name != NULL){
        nopt++;
    }

    struct option * longopts = calloc(nopt+1, sizeof(struct option));
    for(i64 kk = 0; kk < nopt; kk++){
        longopts[kk].name = options[kk].name;
        longopts[kk].has_arg = options[kk].has_arg;
        longopts[kk].flag = options[kk].flag;
        longopts[kk].val = options[kk].val;
    }

    int ch;
    while((ch = getopt_long(argc, argv,
                            "br:c:hmMp:g:v:V",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'b':
            s->write_backbone = 1;
            break;
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
            printf("%s", con2mflock_changelog_txt);
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

    while( L[p] != '\n' && L[p] != '\0'){
        if(L[p] == ','){
            pos[ipos++] = p+1;
            L[p] = '\0';
        }
        p++;
    }

    if(ipos != 5){
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
    if(L[0] == '#'){
        return 1;
    }

    i64 p = 0;
    i64 ipos = 1;
    i64 pos[8];
    pos[0] = 0;

    while( L[p] != '\n' && L[p] != '\0'){
        if(L[p] == '\t'){
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

typedef enum {
    parse_chr_header,
    parse_chr_size
} parse_chr_state;

static void
trim_whitespace(char * str)
{
    if(str == NULL){
        return;
    }
    size_t n = strlen(str);
    size_t first = 0;
    size_t last = n;
    if(n == 0){
        return;
    }

    // Look for whitespaces from the beginning
    char * start = str;
    while(*start != '\0'){
        if(*start == ' '){
            first++;
            start++;
        } else {
            break;
        }
    }

    // Look for whitespaces from the end
    char * end = str+n-1;
    while(end >= start){
        char token = end[0];
        int ws = 0;
        switch(token){
        case ' ':
            ws = 1;
            break;
        case '\n':
            ws = 1;
            break;
        case '\r':
            ws = 1;
            break;
        }
        if(ws){
            last--;
            end--;
        } else {
            break;
        }
    }

    // copy from first to last
    //printf("first=%zu last=%zu\n", first, last);
    size_t wpos = 0;
    for(size_t kk = first; kk < last; kk++){
        str[wpos++] = str[kk];
    }
    str[wpos] = '\0';
    return;
}


static int
parse_chr_definition_header(const char * _line,
                            int * chr_name_col,
                            int * chr_size_col)
{
    char * line = strdup(_line);
    int col = 0;
    char * sub = NULL;
    int has_chr_name = 0;
    int has_chr_size = 0;
    while( (sub = strsep(&line, ",")) ){
        trim_whitespace(sub);
        if(strcmp(sub, "chr_name") == 0){
            *chr_name_col = col;
            has_chr_name = 1;
        }
        if(strcmp(sub, "chr_size") == 0){
            *chr_size_col = col;
            has_chr_size = 1;
        }
        col++;
    }
    free(line);

    if(has_chr_size*has_chr_name == 1){
        return 0;
    }
    return -1;
}

static int
parse_chr_name(const char * sub)
{
    // Excepts an 'r' before the number of letter
    return parse_chr_id(sub);
}

static int
parse_chr_definition_data(const char * _line,
                          int chr_name_col,
                          int * chr_id,
                          int chr_size_col,
                          int * chr_size)
{
    //printf("chr_name_col: %d, chr_size_col: %d\n", chr_name_col, chr_size_col);
    char * line = strdup(_line);
    int col = 0;
    char * sub = NULL;
    int got_chr = 0;
    int got_size = 0;
    while( (sub = strsep(&line, ",")) ){
        //printf("sub%d: %s\n", col, sub);
        if(col == chr_name_col){
            *chr_id = parse_chr_name(sub);
            got_chr = 1;
        }
        if(col == chr_size_col){
            *chr_size = atoi(sub);
            got_size = 1;
        }
        col++;
    }
    free(line);

    if(got_chr*got_size == 1) {
        return 0;
    } else {
        printf("Unparsable line: '%s'\n", _line);
        return -1;
    }
}

// Parse the s->genome file (csv)
// which needs to have the column "chr_name", and "chr_size"
static int get_chr_size(opts * s)
{
    if(s->genome == NULL){
        fprintf(stderr, "--genome not set\n");
        exit(EXIT_FAILURE);
    }

    FILE * fid = fopen(s->genome, "r");

    if(fid == NULL){
        fprintf(stderr, "failed to open %s\n", s->genome);
        perror("get_chr_size");
        exit(EXIT_FAILURE);
    }

    char * line = NULL;
    size_t len = 0;
    ssize_t nread;

    int chr_name_col = -1;
    int chr_size_col = -1;
    int chr_id, chr_size;

    parse_chr_state state = parse_chr_header;
    while ((nread = getline(&line, &len, fid)) != -1) {
        switch(state){
        case parse_chr_header:
            if (nread < 3){
                continue;
            }
            if(line[0] == '#'){
                continue;
            }
            if(parse_chr_definition_header(line, &chr_name_col, &chr_size_col)) {
                printf("%s is not a valid csv header\n", line);;
                goto quit1;
            } else {
                state = parse_chr_size;
            }
            break;
        case parse_chr_size:
            if(nread < 3)
            {
                continue;
            }
            if(line[0] == '#'){
                continue;
            }
            if(parse_chr_definition_data(line,
                                         chr_name_col, &chr_id,
                                         chr_size_col, &chr_size)) {
                printf("Unable to parse '%s'\n", line);
            } else {
                if(s->verbose > 0) {
                    printf("chr_id = %d, chr_size = %d\n", chr_id, chr_size);
                }
                s->chr_size_bp[chr_id] = chr_size;
                if(1 + chr_id > s->n_chr)
                {
                    s->n_chr = 1+chr_id;
                }
            }
            break;
        }
    }

    if(s->verbose > 0){
        printf("s->n_chr = %ld\n", s->n_chr);
        for(int kk = 0; kk < s->n_chr; kk++){
            printf("Chr #%d : %ld\n", kk, s->chr_size_bp[kk]);
        }
    }

    return 0;

quit1:
    fclose(fid);
    free(line);
    return -1;
}


typedef struct {
    char * fname;
    ddict * chr_size;
    u64 n_reads;
    u64 n_lines; // number of lines in the file
} pairs_file_info;

void pairs_file_info_free(pairs_file_info * pinfo){
    if(pinfo == NULL){
        return;
    }
    free(pinfo->fname);
    ddict_free(pinfo->chr_size);
    free(pinfo);
    return;
}

// If the line has the format '#chromosome: chr2 242696752' then write
// 'chr2' to *chr_name and set chr_size to 242696752 and return 1
// If it does not seem like such a line, return 0
//
// On success (return 1) the caller is responsible for freeing the chromosome name
int pfile_get_chromosome(const char * line, char ** _chr_name, u64 * _chr_size)
{
    if(strncmp(line, "#chromosome: ", 13) != 0){
        return 0;
    }
    char * chr_name = strdup(line+13);
    size_t idx = 0;
    while(chr_name[idx] != '\0' && chr_name[idx] != ' '){
        idx++;
    }
    chr_name[idx] = '\0';
    char * chr_size = chr_name + idx + 1;
    *_chr_name = chr_name;
    *_chr_size = atol(chr_size);
    // printf("chr_name = '%s', chr_pos = %ld\n", chr_name, atol(chr_pos));
    return 1;
}


void init_bins_and_labels(opts * s, pairs_file_info * pinfo)
{
    #if 0
    // Set the size for each chromosme
    // and get the total number of bins
    for(i64 kk = 0; kk < s->n_chr; kk++){
        s->chr_size_bin[kk] = (s->chr_size_bp[kk] + s->resolution - 1) / s->resolution;
        s->n_bin += s->chr_size_bin[kk];
    }
#else
    i64 n_chr = ddict_size(pinfo->chr_size);;
    for(i64 kk = 0; kk < n_chr; kk++){
        i64 chr_size = (i64) pinfo->chr_size->entries[kk+1].value;
        s->chr_size_bp[kk] = chr_size;
        s->chr_size_bin[kk] = (chr_size + s->resolution - 1) / s->resolution;
        s->n_bin += s->chr_size_bin[kk];
    }
#endif
    // Create the dense label array for chromflock
    s->labels = malloc(s->n_bin*sizeof(u8));
    size_t idx = 0;
    for(i64 cc = 0; cc < s->n_chr; cc++){
        u8 chr = cc+1;
        i64 n = s->chr_size_bin[cc];
        for(i64 nn = 0 ; nn < n; nn++){
            s->labels[idx++] = chr;
        }
    }
}

// Expecting a pairs1.0 file. Can be gzip'ed
// return the number of lines in s->infile
// and sets s->chr_reads
// uses parse_chr_id (expectations on the chromosome names)
static pairs_file_info *
pairs_file_info_new(const char * fname, int verbose)
{
    if(verbose > 0){
        printf("pairs_info_from_file(%s)\n", fname);
    }

    gzl_state * gzl = gzl_open(fname, 1024);

    if(gzl == NULL){
        printf("pairs_info_from_file: Failed to read %s\n", fname);
        return NULL;
    }

    pairs_file_info * pinfo = calloc(1, sizeof(pairs_file_info));
    pinfo->fname = strdup(fname);
    pinfo->chr_size = ddict_new(1);

    char * L = NULL;
    i64 line_number = 0;
    int gzl_error = 0;

    while( (L = (char *) gzl_get_line(gzl, &gzl_error) ) ){
        line_number++;
        i64 chr1, chr2, pos1, pos2;

        if(L[0] == '#'){
            char * chr_name;
            u64 chr_size;
            // try to parse '#chromosome: chr1 248387328' and similar lines
            if(pfile_get_chromosome(L, &chr_name, &chr_size)){
                ddict_add(pinfo->chr_size, chr_name, (void*) chr_size);
            }
            continue;
        } else {
            if(parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2)){
                printf("pairs_file_info: Unable to parse a line #%ld\n", line_number);
                continue;
            } else {
                pinfo->n_reads++;
            }
        }
    }
    pinfo->n_lines = line_number;
    gzl_destroy(gzl);
    free(L);
    return pinfo;
}

// Write the label array as uint8
// chromflock is happy if we leave 0 unused
static void
write_labels(const opts * s)
{
    if(s->verbose > 0) {
        printf("Writing labels to %s\n", s->label_file);
    }

    if(write_bead_labels(s->label_file, s->labels, s->n_bin)){
        printf("Failed to write to %s\n", s->label_file);
        exit(EXIT_FAILURE);
    }

    return;
}

int cmp_u32_pair(const void * _A, const void * _B)
{
    u32 * A = (u32*) _A;
    u32 * B = (u32*) _B;

    if(A[0] < B[0]){
        return -1;
    }

    if(A[0] > B[0]){
        return 1;
    }

    if(A[1] < B[1]){
        return -1;
    }

    if(A[1] > B[1]){
        return 1;
    }

    return 0;
}

int cmp_rcontact(const void * _A, const void * _B)
{
    rcontact * A = (rcontact*) _A;
    rcontact * B = (rcontact*) _B;

    // first key: chr id
    if(A->chrA < B->chrA){
        return -1;
    }

    if(A->chrA > B->chrA){
        return 1;
    }
    // A->chrA == B->chrA

    // 2nd key: position
    if(A->posA < B->posA){
        return -1;
    }

    if(A->posA > B->posA){
        return 1;
    }

    // A->chrA == B->chrA
    // A->posA == B->posA

    // 3rd key: chrB
    if(A->chrB < B->chrB){
        return -1;
    }

    if(A->chrB > B->chrB){
        return 1;
    }

    // 4th key: posB
    if(A->posB < B->posB){
        return -1;
    }

    if(A->posB > B->posB){
        return 1;
    }

    return 0;
}


static void
write_contacts(opts * s)
{
    assert(s->n_chr > 0);
    assert(s->nlines > 0);

    if(s->verbose > 0){
        printf("Writing contacts to %s\n", s->contact_file);
    }

    gzl_state * gzl = gzl_open(s->infile, 1024);

    if(gzl == NULL){
        gzl_destroy(gzl);
        printf("Failed to read %s\n", s->infile);
        exit(EXIT_FAILURE);
    }

    s->contacts = NULL;
    size_t n_cont_alloc = n_cont_alloc = s->n_reads;
    if(s->write_backbone){
        n_cont_alloc += s->n_bin; // overshoots
    }
    s->contacts = calloc(2*n_cont_alloc, sizeof(u32));

    size_t nchr = s->n_chr;
    if(nchr >= 50){
        exit(EXIT_FAILURE);
    }
    i64 * bin_start = calloc(nchr+1, sizeof(i64));

    for(i64 kk = 1; kk < s->n_chr; kk++){
        bin_start[kk] = bin_start[kk-1] + s->chr_size_bin[kk-1];
    }

    if(s->verbose > 1){
        for(i64 kk = 0; kk < s->n_chr; kk++){
            printf("chr #%ld starts at %ld\n", kk+1, bin_start[kk]);
        }
    }

    i64 line = 0;
    i64 contact_id = 0;
    char * L = NULL;
    int gzl_error;

    while( (L = (char*) gzl_get_line(gzl, &gzl_error) ) ){
        line++;
        i64 chr1, chr2, pos1, pos2;

        if(L[0] == '#'){
            continue;
        }
        parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2);

        if((chr1 > s->n_chr) || (chr2 > s->n_chr)){
            printf("Input error on line %ld\n", line);
            printf("chr1=%ld, chr2=%ld\n", chr1, chr2);
            assert(0);
        }

        i64 bead1 = bin_start[chr1] + pos1/s->resolution;
        i64 bead2 = bin_start[chr2] + pos2/s->resolution;

        if(contact_id >= n_cont_alloc){
            printf("Error: Got more contacts that expected (%ld)\n", s->n_reads);
            printf("While reading line %ld\n", line);
            exit(EXIT_FAILURE);
        }
        s->contacts[2*contact_id + 0] = bead1;
        s->contacts[2*contact_id + 1] = bead2;
        contact_id++;
    }
    gzl_destroy(gzl);


    // Add backbone contacts connecting adjacent beads within the
    // same chromosome if --backbone was provided
    if(s->write_backbone){
        for(i64 kk = 0; kk+1 < s->n_bin; kk++){
            if(s->labels[kk] == s->labels[kk+1]){
                contact_id++;
                s->contacts[2*contact_id + 0] = kk;
                s->contacts[2*contact_id + 1] = kk+1;
            }
        }
    }
    if(s->verbose > 0) {
        printf("Sorting contacts\n");
    }
    qsort(s->contacts, s->n_reads, 2*sizeof(u32), cmp_u32_pair);

    if(s->verbose > 0){
        printf("Removing duplicates\n");
    }
    i64 wpos = 0;
    for(i64 kk = 1; kk < s->n_reads; kk++){
        if(cmp_u32_pair(s->contacts + 2*kk, s->contacts + 2*kk-2) != 0){
            memcpy(s->contacts+2*wpos, s->contacts+2*kk, 2*sizeof(u32));
            wpos++;
        }
    }
    if(s->verbose > 0) {
        printf("Keeping %ld / %ld\n", wpos, s->n_reads);
        printf("Writing contacts to %s\n", s->contact_file);
    }

    if(write_u32(s->contact_file, s->contacts, 2, wpos)){
        printf("Could not write to %s\n", s->contact_file);
        exit(EXIT_FAILURE);
    }

    free(L);
    free(bin_start);
    return;
}

// Write "raw" contacts in the format [[chr_id, pos, chr_id, pos],
// ... ] i.e. without referring to any specific binning
static void
write_contacts_raw(opts * s)
{
    assert(s->n_chr > 0);
    assert(s->nlines > 0);

    if(s->verbose > 0){
        printf("Writing contacts to %s\n", s->contact_file);
    }

    gzl_state * gzl = gzl_open(s->infile, 1024);

    if(gzl == NULL){
        gzl_destroy(gzl);
        printf("Failed to read %s\n", s->infile);
        exit(EXIT_FAILURE);
    }

    size_t n_cont_alloc = n_cont_alloc = s->n_reads*2;
    rcontact * contacts = calloc(n_cont_alloc, sizeof(rcontact));

    i64 line = 0;
    i64 contact_id = 0;
    char * L = NULL;
    int gzl_error;

    while( (L = (char*) gzl_get_line(gzl, &gzl_error) ) ){
        line++;
        i64 chr1, chr2, pos1, pos2;
        if(s->pairs_format){
            if(L[0] == '#'){
                continue;
            }
            parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2);
        } else {
            if(parse_con_line(L, &chr1, &pos1, &chr2, &pos2)){
                printf("Unable to parse a line #%ld\n", line+1);
                exit(EXIT_FAILURE);
            }
        }
        if((chr1 > s->n_chr) || (chr2 > s->n_chr)){
            printf("Input error on line %ld\n", line);
            printf("chr1=%ld, chr2=%ld\n", chr1, chr2);
            assert(0);
        }

        contacts[contact_id].chrA = chr1;
        contacts[contact_id].posA = pos1;
        contacts[contact_id].chrB = chr2;
        contacts[contact_id].posB = pos2;
        contact_id++;
    }
    gzl_destroy(gzl);

    const i64 n_contact = contact_id;
    if(s->verbose > 0) {
        printf("Sorting contacts\n");
    }
    qsort(contacts, n_contact, sizeof(rcontact), cmp_rcontact);

    if(s->verbose > 0){
        printf("Removing duplicates\n");
    }
    i64 wpos = 0;
    for(i64 kk = 1; kk < n_contact; kk++){
        if(cmp_rcontact(contacts + kk-1, contacts + kk) != 0){
            memcpy(s->contacts+wpos, s->contacts+kk, sizeof(rcontact));
            wpos++;
        }
    }

    const i64 n_unique_contact = wpos;
    if(s->verbose > 0) {
        printf("Keeping %ld / %ld\n", n_unique_contact, n_contact);
        printf("Writing contacts to %s\n", s->contact_file_raw);
    }

    if(write_u32(s->contact_file_raw, (u32*) contacts, 4, wpos)){
        printf("Could not write to %s\n", s->contact_file_raw);
        exit(EXIT_FAILURE);
    }

    free(contacts);
    return;
}


static void write_matrix(opts * s)
{
    if(s->verbose > 0){
        printf("Writing contact map to %s\n", s->matrix_file);
    }

    size_t nchr = s->n_chr;
    if(nchr >= 50){
        exit(EXIT_FAILURE);
    }
    i64 * bin_start = calloc(nchr, sizeof(i64));

    for(i64 kk = 1; kk < s->n_chr; kk++){
        bin_start[kk] = bin_start[kk-1] + s->chr_size_bin[kk-1];
    }

    i64 N = s->n_bin;

    if(s->verbose > 0) {
        printf("Matrix size: %ld x %ld\n", N, N);
    }

    u32 * M = calloc(N*N, sizeof(u32));

    gzl_state * gzl = gzl_open(s->infile, 1024);

    if(gzl == NULL){
        printf("Failed to read %s\n", s->infile);
        exit(EXIT_FAILURE);
    }

    i64 line = 0;
    char * L = NULL;
    int gzl_error;
    while( (L = (char *) gzl_get_line(gzl, &gzl_error) ) ){
        line++;
        i64 chr1, chr2, pos1, pos2;
        if(s->pairs_format){
            if(L[0] == '#'){
                continue;
            }
            parse_pairs_line(L, &chr1, &pos1, &chr2, &pos2);
        } else {
            if(parse_con_line(L, &chr1, &pos1, &chr2, &pos2)){
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

    for(i64 kk = 0; kk < N; kk++){
        for(i64 ll = kk +1; ll < N; ll++){
            i64 reads = M[kk + N*ll] + M[ll + N*kk];
            M[kk + N*ll] = reads;
            M[ll + N*kk] = reads;
        }
    }
    if(s->verbose > 0) {
        printf("Writing to %s\n", s->matrix_file);
    }
    if(write_u32(s->matrix_file, M, N, N)){
        printf("Failed to write to %s!\n", s->matrix_file);
        exit(EXIT_FAILURE);
    }

    free(M);
    free(L);
    free(bin_start);
    return;
}

void pairs_file_info_show(pairs_file_info * pinfo, FILE * fid){

    fprintf(fid, "pairs_file_info_show\n");
    fprintf(fid, "   File: %s\n", pinfo->fname);
    int nchr = ddict_size(pinfo->chr_size);
    fprintf(fid, "   Number of lines %ld\n", pinfo->n_lines);
    fprintf(fid, "   Number of reads %ld\n", pinfo->n_reads);
    fprintf(fid, "   Number of chromosomes: %d\n", nchr);

    for(int kk = 1; kk <= nchr; kk++){
        ddict_entry * ent = pinfo->chr_size->entries + kk;
        fprintf(fid, "   %d: %s, %ld\n", kk, (char*) ent->key, (i64) ent->value);
    }
}

void pairs_file_info_write_genome_file(pairs_file_info * pinfo, const char * fname)
{
    FILE * fid = fopen(fname, "w");
    if(fid == NULL){
        printf("Error: Unable to open %s for writing\n", fname);
        return;
    }

    fprintf(fid, "# File generated by chromflock\n");
    fprintf(fid, "# Input data: %s\n", pinfo->fname);
    fprintf(fid, "\n");
    fprintf(fid, "chr_name, chr_size\n");
    int nchr = ddict_size(pinfo->chr_size);
    for(int kk = 1; kk <= nchr; kk++){
        ddict_entry * ent = pinfo->chr_size->entries + kk;
        fprintf(fid, "%s, %ld\n", (char*) ent->key, (i64) ent->value);
    }

    fclose(fid);
}

void validate_chr_mapping(opts * s, pairs_file_info * pinfo)
{
    i64 n_chr = ddict_size(pinfo->chr_size);;
    for(i64 kk = 0; kk < n_chr; kk++){
        const char * chr_name = pinfo->chr_size->entries[kk+1].key;
        i64 chr_id = parse_chr_id(chr_name);
        if(kk != chr_id){
            printf("The pairs file ordering of the chromosomes does not agree with"
                   "how chromosome names are converted to numeric values."
                   "This is simple to fix with a dictionary but wasn't needed until now...\n");;
            printf("%ld != %ld (%s)\n", kk, chr_id, chr_name);

        }
    }
}

int con2mflock(int argc, char ** argv)
{
    opts * s = opts_new();

    parse_command_line(argc, argv, s);

#if 1 // make the program use the dict and this can be skipped....
    // TODO: This should be optional and just override what is in the pairs file
    if(get_chr_size(s)){
        printf("Unable to determine the chromosome sizes from %s\n",
               s->genome);
        exit(EXIT_FAILURE);
    }
#endif

    // Check how many lines of input
    // and how many reads per chromosome
    pairs_file_info * pinfo = pairs_file_info_new(s->infile, s->verbose);
    s->n_reads = pinfo->n_reads;

    if(s->verbose > 0){
        pairs_file_info_show(pinfo, stdout);
    }

    pairs_file_info_write_genome_file(pinfo, s->genome_out_file);

    validate_chr_mapping(s, pinfo);

    // TODO: Agree on the chromosome sizes. I.e.
    // set s->chr_size_bp and what chromosomes to use ...

    // Uses s->resolution and s->chr_size_bp
    // sets
    // s->labels, s->chr_size_bin s->n_bin
    init_bins_and_labels(s, pinfo);

    // with fixed chr sizes we have not even touched the input
    // data so while this works write_contacts might fail leaving us
    // with an incomplete set of output files.
    write_labels(s);

    // - parses the input data for contacts
    // - sorts the contacts
    // - removes duplicates (the lower the resolution, the more duplicates)
    // - writes in a format that chromflock can read
    write_contacts(s);

    // Simply read the contacts, sort them, remove duplicates and
    // write as npy file
    write_contacts_raw(s);

    // Generate a dense contact matrix, not advised for high
    // resolutions
    if(s->gen_matrix){
        write_matrix(s);
    }

    if(s->verbose > 0) {
        printf("done\n");
    }
    opts_free(s);
    pairs_file_info_free(pinfo);
    return EXIT_SUCCESS;
}
