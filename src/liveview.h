#ifndef __view_h_
#define __view_h_

/* Example:
 
#include <liveview.h>
...

  liveXLview xlview;
  xlview.X = X;
  xlview.L = L;
  xlview.N = N;
  xlview.quit = 0;
  xlview.r0 = 0.04; // bead radius

  pthread_t th;
  pthread_create(&th, // thread
      NULL, // pthread_attrib_t
      liveview_t, // function
      &xlview); // arg

// Do some stuff that updates X
//
  xlview.quit = 1;
  pthread_join(th, NULL);
*/


typedef struct {
  size_t N; // Number of points
  double * X; // Coordinates, 3XN
  uint8_t * L; // N Labels
  double r0; // bead radius
  int quit; // Closes when quit = 1;
  elli * E;
} liveXLview;

void * liveview_t(void * conf);


// Static version
int liveview(double * X, uint8_t * L, size_t N, int * quit, double radius, elli * E);

static void rot_x(double *, double);
static void rot_y(double *, double);
static void rot_z(double *, double);
static void matmul(double * , double *, double *);
// static void matshow(double * );

#endif
