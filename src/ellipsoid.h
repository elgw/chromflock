#pragma once

/* Utility functions to handle ellipsoidal geometry
 *
 * References:
 * Feltens J (2009) Vector method to compute the Cartesian (X, Y , Z)
 * to geodetic (, λ, h) transformation on a triaxial ellipsoid. J
 * Geod 83:129–137 10.1007/s00190-008-0246-5
 *
 * Sebahattin Bektas (2014), Orthogonal distance from an ellipsoid,
 * http://dx.doi.org/10.1590/S1982-21702014000400053
 */


#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct{
  double a;
  double b;
  double c;
  double a2;
  double b2;
  double c2;
  size_t maxiter; // Max number of Newton Raphson iterations for gdist
  double eps;
} elli;


/* Create a new ellipsoid. Can be free'd with free */
elli * elli_new(double a, double b, double c);

/* Print it */
void elli_show(elli * E);

/*  */
double elli_dist(elli * E, double * P);

/* Volume */
double elli_vol(elli * E);

/* 1 if point p is inside E, 0 otherwise */
int elli_isInside(elli * E, double * p);

/* Distance to surface */
double elli_dist(elli * E, double * p);

/* Generalized radius
 * if pe is the point where the ray from 0, through p intersects
 * the ellipsoid, ||p||/||pe|| is returned
 */
double elli_radius(elli * E, double * p);

/* Geodesic (=shortest) distance from P to ellipse
 * if Y != NULL, it will return the closest point */
double elli_gdist(elli * ellipse, double * P, double * Y);

/* Using Lagrange multiplier approach */
double elli_gdistL(elli * ellipse,
                   const double * P,
                   double * Y);

/* get scaling parameter but do not appy it to X */
double elli_getScale(const elli * restrict E,
                     const double * restrict X);

/* Same as above but squared (i.e. not sqrt -- faster) */
double elli_getScale2(const elli * restrict E,
                      const double * restrict X);

/* Scale X so that it is on E */
void elli_scale(elli * E, double * X);

/* Calculate the normal at P.
 * Undefined behaviour if P isn't on E
 */
void elli_normal(elli * E, double * P, double * N);

/* Normalize N to 1 */
void elli_nnormal(elli * E, double * P, double * N);
