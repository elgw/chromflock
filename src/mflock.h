#ifndef __mflock_h__
#define __mflock_h__

#include <stdint.h>

typedef struct {
  uint8_t * W; // Matrix with contact indications
  uint32_t * I; // List with pairwise distances
  size_t NI; // Number of pairs in I
  size_t N; // number of points

  uint8_t * L; // chr labels per bead

  double kVol; // Volume exclusion/repulsion of beads
  double kSph; // Spherical confinement
  double kInt; // Interaction of pairs
  double kRad; // Radial force (GPSeq)
  double r0; // bead radius (sphere radius is 1)

  double * R; // wanted radii together with kRad

  // Initialization
  size_t rseed;

  // For the optimizer
  double errstop;
  double gradstop;
  double linestop;
  size_t maxiter;
  size_t maxtime;

  // From the optimizer
  double grad_final;
  size_t iter_final;
  size_t time_final;
  double err_final;

  // Input file names
  char * wfname; // File with contact indications
  char * rfname; // Wanted radius (GPSeq)
  char * xfname; // Initial X
  char * lfname; // Chromosome labels

  // Output file names
  char * xoutfname;
  char * ofoldername; // Outfolder
  char * logfname; // log file name
  int cmmz;

  FILE * logf;

  // general
  int verbose;
  int dynamic;
  int compress; // compress chromosomes by attracting them to their COMs
  int liveView;

} optparam;

void param_show(optparam * p, FILE * f);

optparam * param_copy(optparam * p);
void param_summary(optparam * p, double * X, FILE * f);
void param_readW(optparam * p);
int param_readR(optparam * p);

/* Upscale the points in X2 by a factor of 2 +/- 1 point */
void param_free(optparam * p);


/* For logging */

// Whenever p->verbose > 0
void logwarn(optparam * p, const char *fmt, ...);

// Whenever p->verbose > 1
void loginfo(optparam * p, const char *fmt, ...);

// Whenever p->verbose > 2
void logdebug(optparam * p, const char *fmt, ...);

#endif
