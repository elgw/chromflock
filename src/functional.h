#pragma once

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
#include "cf_util.h"

typedef uint32_t u32;
typedef double f64;


// TODO: Re-use buffers!
typedef struct {
    int nDiv; // number of buckets per dimension
    size_t nH; // number of buckets [nDiv^3]
    u32 * S; // elements per bucket [nH]
    u32 * B; // buckets [nH]
    u32 * C; // start position of each bucket [nH+1]
    f64 * E; // oh, we actually move the points, and not just references to them
    // actually makes sense, according to benchmarks.
} countsort_buffers;

// Settings for the energy landscape
typedef struct {
    double r0; /**< bead radius */
    double dInteraction; /**<  Wanted interaction distance */
    double kVol; /**< Force prohibiting volumetric overlap */
    double kDom; /**< The force that keeps the beads in the domain */
    double kInt; /**< Attraction force for bead pairs */
    double kRad; /**< For radial constraints, using R */
    double kBeadWell; /**< Force that attract specific beads to specific locations */
    double kChrWell; /**< Force that attract beads in specific chromosomes to specific locations */

    mflock_geometry_type geometry;

    elli * E; /**< Ellipsoidal geometry if non-NULL */
    elli * Es; /**<  smaller ellipse defining a "safe region" */

    size_t nIPairs; /**< Number of interaction pairs */
    // TODO: add temporary buffers here
    // TODO: add pointers to data here
    // TODO: add nX as well (which should be called nbead )

    // Position of top and bottom planes that restricts the domain.
    // When set to >1, respectively <-1 these are ignored. The strenght
    // will be kDom
    double top_plane;
    double bottom_plane;
    int diploid;
} mflock_func_t;

/**
 * @breif The total "error" or energy for a certain configuration
 * determined by X
 *
 * @param X the 3D coordinates of the beads
 * @param nX number of beads
 * @param R preferred radial positions
 * @param P list of connected beads (pairs)
 * @param C the settings
 */
double err3(
    const double * restrict X,
    const size_t nX,
    const double * restrict R,
    const uint32_t * restrict P, // list of pairs of interaction
    const mflock_func_t * restrict C
    );

/**
 * @brief The gradient of err3
 *
 * @param X the 3D coordinates of the beads
 * @param nX number of beads
 * @param R preferred radial positions
 * @param P list of connected beads (pairs)
 * @param C the settings
 * @param[out] G Will be set to the gradient with nX elements.
 */
void grad3(
    const double * restrict X,
    const size_t nX,
    const double * restrict R,
    const uint32_t * restrict I, // List of interactions pairs
    double * restrict G,
    const mflock_func_t * restrict C
    );

/** @brief Attract beads to designated wells
 *
 * The wells have a Gaussian shape
 * W(r) = 1 - G(r, sigma)
 * where G is a Gaussian shape scaled so G(0) = 1;
 * r is the distance between the well and a bead
 * sigma is set to fconf->r0
 *
 * This function calculates - \partial W(r) \ partial x_i
 *
 * The wells should be small enough that they don't capture more than
 * one bead. And big enough to increase the capture probability.
 *
 * @param W wells [4 x nW]
 * @param nW
 * @param[out] G gradient vector, will be written to
 */
void
bead_wells_gradient(const mflock_func_t * restrict fconf,
                    const size_t n_bead,
                    const wpos * restrict W,
                    const size_t nW,
                    const double * restrict X,
                    double * restrict G);

/** @brief Attract beads to designated wells
 *
 * The wells have a Gaussian shape
 * W(r) = 1 - G(r, sigma)
 * where G is s Gaussian normalized so G(0) = 1.
 * r is the distance between the well and a bead
 * sigma is set to be fconf->r0
 * W(0) = 0; W(inf) = 1;
 *
 * @param W wells [4 x nW]
 * @param nW
 * @return sum of W(r) for the wells.
 */
double bead_wells_error(const mflock_func_t * restrict fconf,
                        const size_t n_bead,
                        const wpos * restrict W,
                        const size_t nW,
                        const double * restrict X);

/**
 * @brief Reference implementation for err3
 */
double err2(
    double * X,
    size_t nX,
    double * R,
    uint32_t * P, // list of pairs of interaction
    mflock_func_t * C
    );

/**
 * @brief Like err2 but using a matrix, A, to define the beads that
 * should be in contact.
 */
double err(
    double * X,
    size_t nX,
    double * R,
    uint8_t * A,
    mflock_func_t * C
    );


/** Older version of grad, kept for reference */
void grad(
    double * X,
    size_t nX,
    double * R,
    uint8_t * A,
    double * G,
    mflock_func_t * C
    );

/** @brief Older version of grad3, kept for reference.
 *
 * Does not know about ellipsoidal geometry.
 */
void grad2(
    double * X,
    size_t nX,
    double * R,
    uint32_t * I, // List of interactions pairs
    double * G,
    mflock_func_t * C
    );


/** Experimental version of grad3 */
void grad4(
    double * restrict X,
    const size_t nX,
    double * restrict R,
    uint32_t * restrict I, // List of interactions pairs
    double * restrict G,
    const mflock_func_t * restrict C
    );

// For testing only, not really public
double errRepulsion(const double * restrict D,
                    const size_t N,
                    const double d);
