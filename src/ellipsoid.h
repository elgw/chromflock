#pragma once

// Utility functions to handle ellipsoidal geometry
//
// References:
//
// Feltens J (2009) Vector method to compute the Cartesian (X, Y, Z)
// to geodetic (, λ, h) transformation on a triaxial ellipsoid. J
// Geod 83:129–137 10.1007/s00190-008-0246-5
//
// Sebahattin Bektas (2014), Orthogonal distance from an ellipsoid,
// http://dx.doi.org/10.1590/S1982-21702014000400053
//
// To read:
// Alexei Yu. Uteshev, Marina V. Goncharova (2018), Point-to-ellipse
// and point-to-ellipsoid distance equation analysis,
// https://doi.org/10.1016/j.cam.2017.07.021

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


// New ellipsoid with default settings
// a, b, c defines the axis length in x, y, z
elli *
elli_new(double a, double b, double c);

// Print it
void
elli_print(FILE * fid, elli * E);

// Volume
double
elli_vol(elli * E);

// Returns 1 if point p is inside E, 0 otherwise
int
elli_isInside(elli * E, double * p);

// Generalized radial position of a point.
// Returns a positive value where 0
// where a value in means the origo, 1 at the surface, and, > 1
// outside of the surface.
//
// Method:
// if pe is the point where the ray from 0, through p intersects
// the ellipsoid, ||p||/||pe|| is returned
double
elli_radius(const elli * E, const double * p);

// The gdist_ functions computes the shortest distance from P to the
// ellipse E. Y will be set to the closest point on the surface
// (ignored if NULL)
//
// For interior points there are potentially infinite
// number of solutions.

// Following Bektas 2015
// Using The Newton-Raphson method.
// See also this library:
// https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf
double
elli_gdist(const elli * ellipse,
           const double * restrict P,
           double * restrict Y);

// Using Lagrange multiplier approach
double
elli_gdistL(const elli * ellipse,
            const double * restrict P,
            double * restrict Y);

// Root finding method, see John C. Hart, 1994
double
elli_gdistR(const elli * ellipse,
            const double * restrict P,
            double * restrict X);


// Using Sampsons approximation
// d = - F(x) / || \nabla F(x) ||
// where F(x) = x^2/a^2 + y^2/b^2 + z^2/c^2 -1
//
// Gives accurate distances when P is close to E
//
// A negative distance is returned if the point is inside
// A positive distance if the point is outside
//
// Ref: Paul D Sampson. Fitting conic sections to “very scattered”
// data: An iterative refinement of the bookstein algorithm.
// Computer graphics and image processing, 1982
//
// See also:
// https://doi.org/10.48550/arXiv.2401.07114
double
elli_gdistS(const elli * ,
            const double * P);

// Projected distance
//
// Projects X = Projection of p onto E
// returns ||X-p||
double
elli_dist_p(const elli * E, const double * restrict p, double * restrict X);

// get scaling parameter but do not appy it to X
// i.e. the value s such that s*X is on E
double
elli_getScale(const elli * restrict E,
              const double * restrict X);

// Same as above but squared (faster, use when squared distance is enough)
double
elli_getScale2(const elli * restrict E,
               const double * restrict X);

// Project X onto E  along the line segment
// starting at (0,0,0) and passing through X
// If X == (0,0,0), it will be left untouched.
void
elli_project(const elli * E, double * X);



// Calculate the (outward pointing) normal of the ellipsoid at the
// point P.  Undefined behavior if P isn't on E
void
elli_normal(const elli * E,
            const double * restrict P,
            double * restrict N);

// Calculate the (outward pointing) normal of the ellipsoid at the
// point P and normalizes it to unit length. Undefined behavior if P
// isn't on E
void
elli_nnormal(const elli * E,
             const double * restrict P,
             double * restrict N);
