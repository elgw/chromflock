#ifndef __mflock_h__
#define __mflock_h__

#include <stdint.h>
#include "ellipsoid.h"
#include "cf_version.h"

typedef struct {
  uint8_t * W; // Matrix with contact indications
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
  char * wfname; // File with contact indications
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

  // If lua is used to set the parameters as a function of the current iteration
  int luaDynamics;
  char * luaDynamicsFile;

} optparam;

void param_show(optparam * p, FILE * f);

optparam * param_copy(optparam * p);
void param_summary(optparam * p, double * X);
void param_readW(optparam * p);
int param_readR(optparam * p);
double * init_X(optparam * p);
int param_readX(optparam * p, double * X);


/* Upscale the points in X2 by a factor of 2 +/- 1 point */
void param_free(optparam * p);


/* For logging */

// General logging function
void logwrite(optparam * p, int level, const char * fmt, ...);

#endif
