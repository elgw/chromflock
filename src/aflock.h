/**
 * @file aflock.h
 * @author Erik Wernersson
 * @date 2020-2023
 */

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
#include "cf_util.h"
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
    //double th_low;
    //double th_high;
    aflock * af;
    cf_structure * flock;
} adstruct;


/* To be passed to each thread (final_tfun)
 * when running in final mode */
typedef struct{
    size_t thread;
    size_t nThreads;
    size_t nStruct;
    size_t nBeads;

    aflock * af;
    cf_structure * flock;
    // TODO: these do not need double, should be uint16_t
    uint16_t * M; /* Counting number of captured contacts for this thread */
    uint16_t * W; /* Sum of all individual W (uint8_t) for this thread */
    double * rprof; /* Sum of radial values for this thread */
} final_tdata;


/** @brief Initialize a new aflock object */
static aflock *
aflock_init(void);

/* deallocate an aflock object */
static void
aflock_free(aflock *);

/** @brief Load the contact probability matrix shared between the structures */
static void
aflock_load_contact_probabilities(aflock * af);

/** @brief Load all structures specified in af->ffname */
static cf_structure *
aflock_load_structures(aflock * af);

/** @brief Sets both the bead radius and the volume quotient based on
 * one of them */
static void aflock_set_bead_radius(aflock * af);

/** @brief Parse the command line arguments
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 */
static int
aflock_parse_command_line(aflock * p, int argc, char ** argv);

/** @brief Show the command line options */
static void
aflock_usage();


/** @brief Load the coordinates of all structures */
static int
aflock_load_coordinates(aflock * af, cf_structure * c);

/** @brief Initialize a cf_structure object */
static void
cf_structure_init(cf_structure * c, size_t nQ, size_t n);

/** @brief free a cf_structure */
static void
cf_structure_free(cf_structure * c);


/** @brief Create one folder per structure

 * Creates a folder per structure and places the contacts-pairs there.
 * This stage only assigns the contacts with theta == 1, hence the
 * contact pairs are identical for all structures.
 */
static void aflock_init_structures(aflock * af, cf_structure * flock);

/** @brief Update structures by assigning more contacts to them
 *
 * Or, possibly shuffle around contacts withing the given threshold range.
 */
static void aflock_update_structures(aflock * af, cf_structure * flock);

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
static void aflock_finalize_structures(aflock * af, cf_structure * flock);

/* @brief Find what beads that might accept a radial contraint
 *
 * Similar to how contacts are handed out. Used when aflock is called
 * with --rpos set.
 *
 * Update/reassign radial constraints among the chromflock structures
 * only radial positions with a probability of being used < th_high, and >= th_low
 * are used
 *
 * Consider all positions in PR where p<th_high and p>= th_low
 * For each beads:
 * hand the constrains to the structures that has the closest radii already
 * similar to how the contact restraints are handled
 *
 * To start with, all coordinates should already be available in flock.
 * Load preferred radial positions R
 * Load probability/proportion of structures having the constraint into PR */
static void flock_updateR(aflock * af, cf_structure * flock, double th_high, double th_low);


/** @brief Assign new contacts and write them to disk
 *
 * Assigns new contact constraints to the structure c based on
 * the activation distance AD and the contact_probabilities in af->A
 *
 * Finally writes the updated contact constraints to disk.
 *
 * */
static void
cf_structure_assign_contacts(cf_structure * c,
                             aflock * af,
                             double * AD);


/* Float comparison for quicksort */
int cmp_float(const void * A, const void * B);


static float eudist3(const float * A, const float * B);
static float norm3(const float * X);



#endif
