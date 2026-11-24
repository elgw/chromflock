#pragma once

// TODO... integrate with con2mflock which is a mess. Could be the base for a
// rewrite of multiscale_cf_vcn.py

#include <stdint.h>
#include <stdio.h>
#include "ddict.h"

typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;

// Chromosome information
typedef struct {
    char * name; // like chr1 or chrX
    u32 size; // in basepairs
} chromosome_t;

// A polymer is an instance of a chromosome. Any location of
// a polymer can be mapped to a bin_id.
// The mapping of all chromosomes to bin_id is surjective onto
//
typedef struct {
    chromosome_t * chr;
    int copy; // 0, 1, ...
    u32 bin_range[2];
    u32 pos_range[2];
    int poly_id; // used when creating labels for chromflock
} polymer_t;

// Here, a genome is just a collection of chromosomes with unique
// names. TODO: Call it chromosome_set? A genome typically defined as
// the complete set of DNA...
typedef struct {
    char * name;
    char * info;

    polymer_t * poly;
    int n_poly;
    u32 n_bin;

    // A dictionary mapping chromosome names to chromosomes
    ddict * cdict; // chromosome name to index in chrs
    chromosome_t ** chrs; //
    int n_chr; // should support > 1260 ...

} genome_t;

// Given a genome and some binning resolution
// we should be able to map a read
// to a bin, (Chr7, 971201) -> 89000 etc.
typedef struct {
    genome_t * genome;
    u32 resolution; // or bin size
    u32 nbin;
} genome_binner_t;

// Given a set of polymers (that map to chromosomes)
// map reads <-> bins
typedef struct {
    genome_t * genome;
} polymer_binner_t;


// Read a tsv file with the
// columns
// chr_name and chr_size
genome_t * genome_from_tsv(const char * fname);
void genome_free(genome_t * G);
void genome_print(const genome_t * G, FILE * fid);

// Set up a genome binner to allow mapping
// reads <-> bins
genome_binner_t * genome_binner_new(const genome_t * G);

// read to bin
i64 genome_binner_r2b(const genome_binner_t *,
                      const char * chr, u64 pos);

// bin to read returns -1 on failure
int genome_binner_b2r(const genome_binner_t *,
                      u64 bin_id,
                      const char ** chr, // out owned by the genome_binner
                      u64 * pos); // out
