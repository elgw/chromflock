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

/**
 * @brief Settings for the energy landscape
*/
typedef struct conf{
    double r0; /**< bead radius */
    double dInteraction; /**<  Wanted interaction distance */
    double kVol; /**< Force prohibiting volumetric overlap */
    double kDom; /**< The force that keeps the beads in the domain */
    double kInt; /**< Attraction force for bead pairs */
    double kRad; /**< For radial constraints, using R */
    double kWell; /**< Force that attract beads to specific locations */
    elli * E; /**< Ellipsoidal geometry if non-NULL */
    elli * Es; /**<  smaller ellipse defining a "safe region" */

    size_t nIPairs; /**< Number of interaction pairs */
    // TODO: add temporary buffers here
    // TODO: add pointers to data here
    // TODO: add nX as well (which should be called nbead )
} conf;

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
    const conf * restrict C
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
    const conf * restrict C
    );


/**
 * @brief Reference implementation for err3
*/
double err2(
    double * X,
    size_t nX,
    double * R,
    uint32_t * P, // list of pairs of interaction
    conf * C
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
    conf * C
    );


/** Older version of grad, kept for reference */
void grad(
    double * X,
    size_t nX,
    double * R,
    uint8_t * A,
    double * G,
    conf * C
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
    conf * C
    );


/** Experimental version of grad3 */
void grad4(
    double * restrict X,
    const size_t nX,
    double * restrict R,
    uint32_t * restrict I, // List of interactions pairs
    double * restrict G,
    const conf * restrict C
    );

// For testing only, not really public
double errRepulsion(const double * restrict D,
                    const size_t N,
                    const double d);

#endif
