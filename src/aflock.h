#ifndef _aflock_h_
#define _aflock_h_

#define _GNU_SOURCE

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
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

#include "contact_pairs_io.h"
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

typedef enum {
    MODE_UNKNOWN, /* Mode not set */
    MODE_UPDATE, /* Update existing structures */
    MODE_INIT, /* First iterations */
    MODE_FINAL, /* Perform some measurements */
    MODE_UPDATE_BLOCKED, /* ?? */
    MODE_EXPERIMENTAL /* ..  */
} aflock_mode;


/* Configuration and state of the program */
typedef struct {
    aflock_mode mode; // What to do

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

    /* The command line arguments to mflock which will be written to mflock-jobs */
    char * mflock_arguments;

    // Flags
    int experimental;
    /* Haploid or diploid */
    bool diploid;

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
    size_t mem_limit;
} aflock;


/* For each individual chromflock structure */
typedef struct {
    char * contact_pairs_file; /* Contact indicator matrix */
    char * xfName; /* Coordinates */
    char * rfName; /* Radial constraints, e.g., GPSeq */

    float * X; /* 3xN 3D coordinates */
    uint32_t * Q; /* Queued contacts to write, 2xQS */
    size_t nQ; /* Number of pairs in the queue */
    uint8_t * W; /* Contact indicator matrix */
} cf_structure;

/* Activation Distance Struct
 * To be passed to threads */
typedef struct{
    double * A; /* For storing activation distances */
    double * AD;
    size_t thread;
    size_t nThreads;
    double th_low;
    double th_high;
    aflock * fc;
    cf_structure * flock;
} adstruct;


/* To be passed to each thread (final_tfun)
 * when running in final mode */
typedef struct{
    size_t thread;
    size_t nThreads;
    size_t nStruct;
    size_t nBeads;

    aflock * fc;
    cf_structure * flock;
    // TODO: these do not need double, should be uint16_t
    double * M; /* Counting number of captured contacts for this thread */
    double * W; /* Sum of all individual W (uint8_t) for this thread */
    double * rprof; /* Sum of radial values for this thread */
} final_tdata;


// Initialize a new aflock object
static aflock * aflock_init(void);
/* deallocate an aflock object */
static void aflock_free(aflock *);

// Load all structures specified in fc->ffname
static cf_structure * load_structures(aflock * fc);
// Initialize a cf_structure object
static void cf_structure_init(cf_structure * c, size_t nQ, size_t n);
// Write a complete W matrix for a structure,
// i.e., no update. Used in MODE_INIT only

static void cf_structure_free(cf_structure * c);
static int cf_structure_load_X(aflock * fc, cf_structure * c);

// The main functionality,

/* Update/reassign radial constraints among the chromflock structures
 * only radial positions with a probability of being used < th_high, and >= th_low
 * are used
 */
static void flock_updateR(aflock * fc, cf_structure * flock, double th_high, double th_low);

/** @brief Assign new contacts and write them to disk
 *
 * Write a new W matrix to disk, setting contacts between beads to
 * 1 if the distance is below the activation distance given by AD
 *
 * This will read an existing W matrix and only update contacts for which
 * A is in the threshold range.
 * */
static void struct_write_W_AD(aflock * fc,
                              cf_structure * c,
                              double * AD,
                              double th_high,
                              double th_low);

// Float comparison for quicksort
int cmp_float(const void * A, const void * B);

/* Load the contact probability matrix shared between the structures */
static void aflock_load_A(aflock * fc);

static float eudist3(const float * A, const float * B);
static float norm3(const float * X);

static int argparsing(aflock * p, int argc, char ** argv);
static void usage();

#endif
