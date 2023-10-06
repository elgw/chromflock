#ifndef _functional_h_
#define _functional_h_

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef MATLAB
#include "mex.h"
#endif

#include "ellipsoid.h"

typedef struct conf{
    double r0;
    double dInteraction; // Wanted interaction distance
    double kVol;
    double kDom;
    double kInt;
    double kRad; // For radial constraints, using R
    size_t nIPairs; // Number of wanted interactions
    elli * E;
    elli * Es; // smaller ellipse defining "safe region"
} conf;

double err(
    double * X,
    size_t nX,
    double * R,
    uint8_t * A,
    conf * C
    );

double err2(
    double * X,
    size_t nX,
    double * R,
    uint32_t * P, // list of pairs of interaction
    conf * C
    );

// With hash structure for collisions
double err3(
    double * restrict X,
    const size_t nX,
    double * restrict R,
    uint32_t * restrict P, // list of pairs of interaction
    const conf * restrict C
    );

void grad(
    double * X,
    size_t nX,
    double * R,
    uint8_t * A,
    double * G,
    conf * C
    );

void grad2(
    double * X,
    size_t nX,
    double * R,
    uint32_t * I, // List of interactions pairs
    double * G,
    conf * C
    );

// hashed only one that works with elliptical geometry
void grad3(
    double * restrict X,
    const size_t nX,
    double * restrict R,
    uint32_t * restrict I, // List of interactions pairs
    double * restrict G,
    const conf * restrict C
    );

// Experimental
void grad4(
    double * restrict X,
    const size_t nX,
    double * restrict R,
    uint32_t * restrict I, // List of interactions pairs
    double * restrict G,
    const conf * restrict C
    );

// For testing only, not really public
double errRepulsion(double * restrict D,
                    const size_t N,
                    const double d);

#endif
