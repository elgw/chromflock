#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <SDL.h>
#include "hsvrgb.h"
#include "ellipsoid.h"

/* Idea:
 * Create a light weight renderer that can be plugged in anywhere to monitor the status of some
 * matrix
 *
 *  - Use pthreads to run separately
 *  - Update window whenever the appointed data is changed
 *    or use some kind of signalling.
 *    Also add a blocking command to destroy the window (before the data is removed).
 *
 *
 * Example:

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
  volatile int quit; // Closes when quit = 1;
  elli * E;
} liveXLview;

void * liveview_t(void * conf);

int
liveview(double * X,
         uint8_t * L,
         size_t N,
         volatile int * quit,
         double radius,
         elli * E);
