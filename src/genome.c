#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "genome.h"

// Forward declarations
static int
genome_add_chr(genome_t * G, const char * chr_name, int chr_size);

void chromosome_free(chromosome_t * C){
    if(C == NULL){
        return;
    }
    free(C->name);
    free(C);
    return;
}

// Removes whitespaces from both ends
// and moves the data so that the pointer does not have to
// be updated. Obviously a bad strategy for large strings.
// Example:
// '\t x\t\n\0' => 'x\0'
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


void genome_free(genome_t * G)
{
    for(u64 kk = 1; kk <= G->cdict->n_entries; kk++){
        chromosome_free(G->cdict->entries[kk].value);
    }
    ddict_free(G->cdict);
    free(G->chrs);
    free(G->name);
    free(G->info);
    free(G);
}

void genome_print(const genome_t * G, FILE * fid){
    assert(G != NULL);
    fprintf(fid, "Genome: %s\n", G->name);
    for(int kk = 0; kk < G->n_chr; kk++){
        fprintf(fid, "%s (%u)\n",
                G->chrs[kk]->name, G->chrs[kk]->size);
    }
    return;
}

typedef enum {
    parse_chr_header,
    parse_chr_size
} parse_chr_state;

static int
parse_chr_definition_header(const char * _line,
                            int * chr_name_col,
                            int * chr_size_col)
{
    char * line = strdup(_line);
    char * line0 = line;
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
    free(line0);

    if(has_chr_size*has_chr_name == 1){
        return 0;
    }
    return -1;
}

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

static int
parse_chr_definition_data(const char * _line,
                          int chr_name_col,
                          char ** chr_name,
                          int chr_size_col,
                          int * chr_size)
{
    //printf("chr_name_col: %d, chr_size_col: %d\n", chr_name_col, chr_size_col);
    char * line = strdup(_line);
    char * line0 = line;
    int col = 0;
    char * sub = NULL;
    int got_chr = 0;
    int got_size = 0;
    while( (sub = strsep(&line, ",")) ){
        //printf("sub%d: %s\n", col, sub);
        if(col == chr_name_col){
            assert(got_chr == 0);
            chr_name[0] = strdup(sub);
            got_chr = 1;
        }
        if(col == chr_size_col){
            *chr_size = atoi(sub);
            got_size = 1;
        }
        col++;
    }
    free(line0);

    if(got_chr*got_size == 1) {
        return 0;
    } else {
        if(got_chr){
            free(chr_name[0]);
        }
        printf("Unparsable line: '%s'\n", _line);
        return -1;
    }
}


// Read lines from fid to populate G
// which is initialized to 0
// return 0 on success
// or 1 on failure, then genome_free(G) is required
static int
parse_genome_file(FILE * fid,
                  genome_t * G,
                  const int verbose)
{
    char * line = NULL;
    size_t len = 0;
    ssize_t nread;

    int chr_name_col = -1;
    int chr_size_col = -1;
    int chr_size = 0;
    char * chr_name = NULL;

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
            if(parse_chr_definition_header(line,
                                           &chr_name_col,
                                           &chr_size_col)) {
                printf("%s is not a valid csv header\n", line);;
                free(line);
                return 1;
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
                                         chr_name_col, &chr_name,
                                         chr_size_col, &chr_size)) {
                printf("Unable to parse '%s'\n", line);
            } else {
                if(verbose > 1) {
                    printf("chr_name = %s, chr_size = %d\n", chr_name, chr_size);
                }

                genome_add_chr(G, chr_name, chr_size);
                free(chr_name);
            }
            break;
        }
    }
    free(line);
    G->n_chr = ddict_size(G->cdict);
    printf("Found %d chromosomes\n", G->n_chr);
    G->chrs = malloc(G->n_chr*sizeof(chromosome_t*));
    for(int kk = 0; kk < G->n_chr; kk++){
        G->chrs[kk] = (chromosome_t *) G->cdict->entries[kk+1].value;
    }

    return 0;
}

// Note: this does not update G->chrs or G->n_chr
static int
genome_add_chr(genome_t * G, const char * chr_name, int chr_size)
{
    if(ddict_get(G->cdict, "chr_name"))
    {
        printf("Error: %s alread in genome\n", chr_name);
        return 1;
    }
    chromosome_t * chr = calloc(1, sizeof(chromosome_t));
    chr->name = strdup(chr_name);
    chr->size = chr_size;
    ddict_add(G->cdict, chr_name, (void *) chr);
    return 0;
}


static genome_t *
genome_from_csv(const char * fname, int verbose)
{
    if(fname == NULL){
        assert(fname != NULL);
        return NULL;
    }

    FILE * fid = fopen(fname, "r");

    if(fid == NULL){
        assert(fid != NULL);
        fprintf(stderr, "failed to open %s\n", fname);
        return NULL;
    }

    genome_t * G = calloc(1, sizeof(genome_t));
    G->cdict = ddict_new(1);

    if(parse_genome_file(fid, G, verbose)){
        fclose(fid);
        genome_free(G);
        return NULL;
    }

    G->name = strdup(fname);

    fclose(fid);
    return G;
}

int main(int argc, char ** argv){

    if(argc < 2){
        printf("Usage: %s genome.csv\n", argv[0]);
        return EXIT_FAILURE;
    }
    int verbose = 1;

    genome_t * G = genome_from_csv(argv[1], verbose);
    genome_print(G, stdout);
    genome_free(G);

    // Next step: Create a binner to be able to parse HiC data. I.e.
    // at a given resolution, map a genomic coordinate like (chr7, 89409)
    // to a bin number.

    genome_t * G = genome_from_csv(argv[1], verbose);
    G = genome_remove_chr(G, "chrY"); // TODO
    int res = 1000000;
    genome_binner_t * GB = genome_binner_new(G, res); // TODO
    printf("At %d bp resolution there are %ld bins\n",
           GB->nbin);
    genome_binner_free(GB); // TODO
    genome_free(G);


    // Next step: Write a table of all the bins

    // Next step: create the collection of polymers -- what is that
    // called? From the genome. For sc_chromflock

    return EXIT_SUCCESS;
}
