#ifndef _aflock_h_
#define _aflock_h_

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

#ifdef __linux__
#include <dlfcn.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#endif

#include "cf_version.h"
#include "ellipsoid.h"
#include "oscp.h"
#include "wio.h"

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

#define MODE_UNKNOWN -1
#define MODE_UPDATE 0
#define MODE_INIT 1
#define MODE_FINAL 2
#define MODE_UPDATE_BLOCKED 3
#define MODE_EXPERIMENTAL 4
#include "ellipsoid.h"

// Meta data etc, flock configuration
typedef struct {
    int mode; // What to do // TODO: enumerate

    char * afname; // File with contact probability matrix
    double * A; // Contact probability matrix

    size_t nBeads; // Number of beads per structure
    size_t nStruct; // Number of structures

    // Use contact within this range from A
    double th_low;
    double th_high;

    // Size of write queue for each structure
    size_t QS;

    // Geometry
    double vq; // volume quotient
    double r0; // bead radius
    double ea;
    double eb;
    double ec;
    elli * E;
    double dContact; // contact threshold

    char * mflock_arguments; // Will be written to mflock-jobs

    // Flags
    int experimental;
    int diploid;

    char * rfname;
    char * prfname;

    /* For DamID
     * R = 1; probR = DamID score;
     * For GPSeq:
     * R = GPSeq radius, probR = 1;, simpler is to use the same R.double for all structures;
     */

    double * R;
    double * probR;
    size_t nThreads;
    int verbose;
} fconf;


/* For each individual chromflock structure */
typedef struct {
    // All points are loaded
    float * X; // 3xN

    uint32_t * Q; // Queued contacts to write, 2xQS
    size_t nQ; // Number of pairs in the queue

    // File names
    char * wfName; // Contact indicator matrix
    char * xfName; // Coordinates
    char * rfName;

    // size_t N;
    uint8_t * W;
} chrom;

// Initialize a new fconf object
static fconf * fconf_init(void);
/* deallocate an fconf object */
static void fconf_free(fconf *);

// Load all structures specified in fc->ffname
static chrom * load_structures(fconf * fc);


// Load all coordinates from all structures
static chrom * load_structures(fconf * fc);

// Initialize a chrom object
static void chrom_init(chrom * c, size_t nQ, size_t n);


// Write a complete W matrix for a structure,
// i.e., no update. Used in MODE_INIT only
static void struct_write_W0(fconf * fc, chrom * cc, uint8_t * W0);

static void ch_free(chrom * c);
static int ch_load_X(fconf * fc, chrom * c);

// The main functionality,

/* Update/reassign radial constraints among the chromatin structures in flock
 * only radial positions with a probability of being used < th_high, and >= th_low
 * are used
 */
static void flock_updateR(fconf * fc, chrom * flock, double th_high, double th_low);


// Float comparison for quicksort
int cmp_float(const void * A, const void * B);

/* Load the contact probability matrix shared between the structures */
static void fconf_load_A(fconf * fc);

static float eudist3(float * A, float * B);
static float norm3(float * X);

static int argparsing(fconf * p, int argc, char ** argv);
static void usage();

#endif
