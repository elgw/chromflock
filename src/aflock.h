#ifndef _aflock_h_
#define _aflock_h_

#define MODE_UNKNOWN -1
#define MODE_UPDATE 0
#define MODE_INIT 1
#define MODE_FINAL 2
#define MODE_UPDATE_BLOCKED 3
#define MODE_EXPERIMENTAL 4
#include "ellipsoid.h"

// Meta data etc, flock configuration
typedef struct {
    int mode; // What to do

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
    int get_quality;
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


// For each individual chromatine structure
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
fconf * fconf_init(void);
// deallocate an fconf object
void fconf_free(fconf *);

// Load all structures specified in fc->ffname
chrom * load_structures(fconf * fc);
void struct_write_W(fconf * fc, chrom * c);
// Load the contact probability matrix shared between the structures
void fconf_load_A(fconf * fc);

// Load all coordinates from all structures
chrom * load_structures(fconf * fc);


// Initialize a chrom object
int chrom_init(chrom * c, size_t nQ, size_t n);
// Assign a contact pair (kk, ll) to a chrom structure
void chrom_assign(fconf * fc, chrom * ch, size_t kk, size_t ll);
// Write down all assigned contact pairs for a structure
// and empty the queue, used in MODE_UPDATE
void struct_write_W(fconf * fc, chrom * c);
// Write a complete W matrix for a structure,
// i.e., no update. Used in MODE_INIT only
void struct_write_W0(fconf * fc, chrom * cc, uint8_t * W0);

void ch_free(chrom * c);
int ch_load_X(fconf * fc, chrom * c);

// The main functionality,
// assign contacts to the most suited structures
void flock_assign(fconf * fc, chrom * flock, size_t kk, size_t ll, float prob);

// Build a contact map from all structures
// for the final mode
void struct_get_interactions(fconf * fc, chrom * cc, double * M, double cDistance);

/* Update/reassign radial constraints among the chromatin structures in flock
 * only radial positions with a probability of being used < th_high, and >= th_low
 * are used
 */
void flock_updateR(fconf * fc, chrom * flock, double th_high, double th_low);


// Float comparison for quicksort
int cmp_float(const void * A, const void * B);


static float eudist3(float * A, float * B);
static float norm3(float * X);

static int argparsing(fconf * p, int argc, char ** argv);
static void usage();

#endif
