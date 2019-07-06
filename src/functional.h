#ifndef _functional_h_
#define _functional_h_

#include <stdint.h>

typedef struct conf{
  double r0;
  double dInteraction; // Wanted interaction distance
  double kVol;
  double kSph;
  double kInt;
  double kRad; // For radial constraints, using R
  size_t nIPairs; // Number of wanted interactions
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

// hashed
void grad3(
    double * restrict X, 
    const size_t nX, 
    double * restrict R, 
    uint32_t * restrict I, // List of interactions pairs
    double * restrict G,
    const conf * restrict C
    );


#endif

