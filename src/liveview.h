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

// SDL2 based bead viewer
//
// Has to be run from the main thread
// Busy waiting/drawing

int
liveview(const double * XYX, // 3xN coordinates -- a pointer into mflock
         const uint8_t * labels, // N labels
         size_t n_bead,
         volatile int * quit, // Set to 1 to quit and deallocate
         double radius, // bead radius
         const elli * E); // set to NULL if spherical domain
