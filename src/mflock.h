/**
 * @file mflock.h
 * @author Erik Wernersson
 * @date 2020-2023
 */

#pragma once

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
#include "fast_prng/normal.h"

#ifdef SDL
#include <pthread.h>
#include "liveview.h"
#include "liveview.c"
#endif

#include "cf_version.h"
#include "cf_util.h"
#include "cmmwrite.h"
#include "ellipsoid.h"
#include "functional.h"
#include "wio.h"
#include "contact_pairs_io.h"

typedef struct {
    uint32_t * I; // List with pairwise distances
    size_t NI; // Number of pairs in I
    size_t N; // number of points

    uint8_t * L; // chr labels per bead
    double * R; // wanted radii together with kRad

    // Geometry
    // Sphere if E isn't set.
    double r0; // bead radius
    double volq;
    elli * E; // Ellipse

    // Forces
    double kVol; // Volume exclusion/repulsion of beads
    double kDom; // Spherical confinement
    double kInt; // Interaction of pairs
    double kRad; // Radial force (GPSeq)


    // Initialization
    size_t rseed;

    // For the optimizer
    size_t maxiter;
    size_t maxtime;

    // From the optimizer
    double grad_final;
    size_t iter_final;
    size_t time_final;
    double err_final;

    // Input file names
    char * wfname; // File with contact indications (depreciated)
    char * contact_pairs_file; /* File with contact pairs */
    char * rfname; // Wanted radius (GPSeq)
    char * xfname; // Initial X
    char * lfname; // Chromosome labels

    int newx;

    // Output file names
    char * xoutfname;
    char * ofoldername; // Outfolder
    char * logfname; // log file name
    int cmmz;

    FILE * logf;

    // general
    int verbose;
    double compress; // compress chromosomes by attracting them to their COMs
    int liveView;

    /* Name of lua script to handle the beads dynamics */
    char * luaDynamicsFile;

} mflock_t;


/** @brief create a new default configuration
 * free with mflock_free
*/
mflock_t *  mflock_new(void);

/** @brief free an mflock_t
 *
 * also frees everything that it points to
 */
void mflock_free(mflock_t * p);

/** @brief Print the settings to FILE
 *
 * f can of course be stdout.
 */
void mflock_show(mflock_t * p, FILE * f);

/** @brief Report status of mflock to log and screen
 *
 */
void mflock_summary(mflock_t * p, const double * restrict X);

/** @brief Read contact pairs from a binary file
    Sets p->I (the contacts) and p->NI (number of contact pairs)
    @return - Nothing, but aborts the program on failure.
*/
static void mflock_read_contact_pairs(mflock_t * p);

/** @brief Load radial constraints
 *
 * only if rfname is set
 * Read GPSeq radius values as binary double.
 * If that does not work, try as text, one value per line
 */
int mflock_load_radial_constraints(mflock_t * p);


double * init_X(mflock_t * p);

/** @brief Load or set new coordinates
 *
 * Tries to call mflock_init_coordinates
 *
*/

double * mflock_init_coordinates(mflock_t * p);


/** @brief Load bead coordinates from csv file
 *
 */
int mflock_load_coordinates(mflock_t * p, double * X);


/** @brief Write coordinates to disk, also write the column names
 * to the log file  */
void mflock_save_coordinates(mflock_t * p, double * X);

/** @brief Read label matrix pointed to by p->lfname
 */
static int mflock_load_bead_labels(mflock_t * p);

/* For logging */
void logwrite(mflock_t * p, int level, const char * fmt, ...);

/**
 * @breif The beads dynamics main loop
 *
 * @param X the coordinates of each bead
 * @param p the settings
 * @param Fb: the Brownian force
 */
static int mflock_dynamics(double * restrict X,
                   mflock_t * restrict p,
                   double Fb);

/** @brief parse command line arguments */
static int mflock_parse_cli(mflock_t * p, int argc, char ** argv);

/** @brief initialization from valid command line arguments
*/
static int mflock_init(mflock_t * p, int argc, char ** argv);


/** @brief double check before running */
static void mflock_validate(mflock_t * p);
