#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>

#ifdef MATLAB
#include "mex.h"
#endif

#include "functional.h"
#define INLINED inline __attribute__((always_inline))

// Consider compiling with -ffinite-math-only

INLINED static double norm3(double * restrict X)
{
  double n = 0;
  for(size_t kk = 0; kk<3; kk++)
    n+=pow(X[kk], 2);
  return sqrt(n);
}

INLINED static double norm32(double * restrict X)
{
  double n = 0;
  for(size_t kk = 0; kk<3; kk++)
    n+=pow(X[kk], 2);
  return n;
}

INLINED static double eudist3sq(double * A, double * B)
{
  /* SQUARED Euclidean distance between two 3D-vectors */
  return pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2],2);
}

INLINED static double eudist3(double * A, double * B)
{
  /* Euclidean distance between two 3D-vectors */
  return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}

INLINED size_t hash_coord(const int nDiv, const double X)
{
  double v = (X+1)/2 * nDiv;

  assert(isfinite(v) == 1);

  if(v<=0)
    return 0;
  if(v>=nDiv)
    return nDiv-1;
  return floor(v);
}


INLINED size_t hash(const size_t nDiv, double * restrict X)
{

  size_t h = hash_coord(nDiv, X[0]) + 
    nDiv*hash_coord(nDiv, X[1]) + 
    nDiv*nDiv*hash_coord(nDiv, X[2]);

  return h;
}


static double errRepulsion(double * restrict D, 
    const size_t N, 
    const double d)
{

  // NOTE: Will crash if D contains nan or inf values.

#ifndef NDEBUG
  // Verify that D was actually allocated
  for(size_t kk = 0; kk<3*N; kk++)
  {
    assert(isfinite(D[kk]));
  }
  fflush(stdout);
#endif


  // Mostly copied from volHASH/src/volhash3.c
  double d2 = pow(d,2);

  // 1. Figure out how many elements per bucket
  int nDiv = 9;
  size_t nH = pow(nDiv, 3);
  uint32_t * S = malloc(nH*sizeof(uint32_t));

  for(size_t kk = 0; kk<nH; kk++)
  {
    S[kk] = 0;
  }

  for(size_t kk = 0; kk<N; kk++)
  {
    S[hash(nDiv, D+3*kk)]++;
  }

  /* Create boundaries for the data */
  uint32_t * B = malloc((nH+1)*sizeof(uint32_t));
  uint32_t * C = malloc((nH+1)*sizeof(uint32_t));

  assert(B != NULL);

  B[0] = 0;
  C[0] = 0; // Start positions when writing
  for(size_t kk = 1; kk<=nH; kk++)
  {
    B[kk] = B[kk-1]+S[kk-1];
    C[kk] = B[kk];
  }

  if(0){
    for(size_t kk = 0; kk<nH; kk++)
    {
      printf("[%u -- %u], S[%zu] = %u\n", B[kk], B[kk+1], kk, S[kk]);
    }
  }

  // Dots sorted according to their bucket
  double * E = malloc(3*N*sizeof(double));

  /* Move data into new structure */
  for(size_t kk = 0; kk<N; kk++)
  {
    size_t h = hash(nDiv, D+3*kk);

    size_t writepos = 3*C[h];
    C[h]++;
    for(int idx = 0; idx<3; idx++)
    {
      assert(writepos<N*3);
      E[writepos++] = D[3*kk+idx];      
    }
  }

  /* 
   * now we can hash(X) a point X,
   * using B to see where the bucket is
   * in E 
   * */

  // Get some job done!
  double err = 0;
  for(size_t kk = 0; kk<N; kk++) 
    /* Looping over E in order to avoid self-matches 
     * and duplicates */
  {
    double deps = d;
    size_t ha_min = hash_coord(nDiv, E[3*kk]-deps);
    size_t ha_max = hash_coord(nDiv, E[3*kk]+deps);

    size_t hb_min = hash_coord(nDiv, E[3*kk+1]-deps);
    size_t hb_max = hash_coord(nDiv, E[3*kk+1]+deps);

    size_t hc_min = hash_coord(nDiv, E[3*kk+2]-deps);
    size_t hc_max = hash_coord(nDiv, E[3*kk+2]+deps);

    for(size_t aa = ha_min; aa <= ha_max; aa++)    
      for(size_t bb = hb_min; bb <= hb_max; bb++)    
        for(size_t cc = hc_min; cc <= hc_max; cc++)
        {

          size_t hash = 
            aa + 
            bb*nDiv + 
            cc*pow(nDiv,2);

          //          printf("%d, %d, %d, hash: %zu / %zu\n", aa, bb, cc, hash, nH); fflush(stdout);

          for(size_t pp = B[hash]; pp<B[hash+1]; pp++)
          {
            if(pp>kk) {
              double dist2 = eudist3sq(E+3*pp, E+3*kk);
              if(dist2<d2)
              {
                err += pow(sqrt(dist2) - d, 2);
              }
            }
          }
        }
  }

  free(B);
  free(C);
  free(E);
  free(S);

  return err;
}


static double gradRepulsion(double * restrict D, 
    double * restrict G, 
    const size_t N, 
    const double d, 
    const double kVol)
{

#ifndef NDEBUG
  // Verify that D was actually allocated
  for(size_t kk = 0; kk<3*N; kk++)
  {
    assert(isfinite(D[kk]));
  }
  fflush(stdout);
#endif

  // Mostly copied from volHASH/src/volhash3.c
  double d2 = pow(d,2);

  // 1. Figure out how many elements per bucket
  int nDiv = 9;
  size_t nH = pow(nDiv, 3);
  uint32_t * S = malloc(nH*sizeof(uint32_t));

  for(size_t kk = 0; kk<nH; kk++)
  {
    S[kk] = 0;
  }

  for(size_t kk = 0; kk<N; kk++)
  {
    S[hash(nDiv, D+3*kk)]++;
  }

  /* Create boundaries for the data */
  uint32_t * B = malloc((nH+1)*sizeof(uint32_t));
  uint32_t * C = malloc((nH+1)*sizeof(uint32_t));

  B[0] = 0;
  C[0] = 0; // Start positions when writing
  for(size_t kk = 1; kk<=nH; kk++)
  {
    B[kk] = B[kk-1]+S[kk-1];
    C[kk] = B[kk];
  }

  if(0){
    for(size_t kk = 0; kk<nH; kk++)
    {
      printf("[%u -- %u], S[%zu] = %u\n", B[kk], B[kk+1], kk, S[kk]);
    }
  }

  // Dots sorted according to their bucket
  double * E = malloc(3*N*sizeof(double));
  size_t * P = malloc(N*sizeof(double)); // Keep also bead numbers

  /* Move data into new structure */
  for(size_t kk = 0; kk<N; kk++)
  {
    size_t h = hash(nDiv, D+3*kk);

    size_t writepos = 3*C[h];
    C[h]++;
    P[writepos/3] = kk;
    for(int idx = 0; idx<3; idx++)
    {
      assert(writepos<N*3);
      E[writepos++] = D[3*kk+idx];      
    }
  }

  /* 
   * now we can hash(X) a point X,
   * using B to see where the bucket is
   * in E 
   * */

  // Get some job done!
  double err = 0;
  for(size_t kk = 0; kk<N; kk++) 
    /* Looping over E in order to avoid self-matches 
     * and duplicates */
  {
    double deps = d;
    size_t ha_min = hash_coord(nDiv, E[3*kk]-deps);
    size_t ha_max = hash_coord(nDiv, E[3*kk]+deps);
    size_t hb_min = hash_coord(nDiv, E[3*kk+1]-deps);
    size_t hb_max = hash_coord(nDiv, E[3*kk+1]+deps);
    size_t hc_min = hash_coord(nDiv, E[3*kk+2]-deps);
    size_t hc_max = hash_coord(nDiv, E[3*kk+2]+deps);

    for(size_t aa = ha_min; aa <= ha_max; aa++)    
      for(size_t bb = hb_min; bb <= hb_max; bb++)    
        for(size_t cc = hc_min; cc <= hc_max; cc++)
        {

          size_t hash = 
            aa + 
            bb*nDiv + 
            cc*pow(nDiv,2);
          for(size_t pp = B[hash]; pp<B[hash+1]; pp++)
          {
            if(pp>kk) {
              double dist2 = eudist3sq(E+3*pp, E+3*kk);
              if(dist2<d2)
              {

                // Retrieve original positions
                size_t KK = P[kk];
                size_t LL = P[pp];

                double di = sqrt(dist2);
                for(int idx = 0; idx<3; idx++)
                {
                  G[3*KK+idx] += kVol*2*(E[3*kk+idx] - E[3*pp+idx])/di*(di - d);
                  G[3*LL+idx] -= kVol*2*(E[3*kk+idx] - E[3*pp+idx])/di*(di - d);
                }

              }
            }
          }
        }
  }


  free(B);
  free(C);
  free(E);
  free(S);
  free(P);

  return err;
}

double err3(double * restrict X, const size_t nX, double * restrict R, uint32_t * restrict P, const conf * restrict C )
{
  /* Alternative version with a list of pairs in contact instead of A */

  // Wanted radii
  double errRad = 0;
  if( (R != NULL) & (C->kRad > 0))
  {
    for(size_t kk = 0; kk<nX; kk++)
    {
      if(isfinite(R[kk]) == 1)
      { // TODO: create list of all finite first, no need to have a branch here
        const double r = norm3(X+kk*3);
        errRad += pow(r-R[kk], 2);
      }
    }
  }

  // Keep inside sphere
  double errSph = 0;
  double ds = (1-C->r0);
  double ds2 = pow(ds, 2);
  for(size_t kk = 0 ; kk<nX; kk++)
  {
    const double r2 = norm32(X+kk*3); // norm^2
    // TODO: use min/max here to avoid branch?
    if( r2 > ds2)
    {
      double r = sqrt(r2);
      errSph += pow(r - ds, 2);
    }
  }

  // Wanted contacts/interactions
  double errInt = 0;
  for(size_t pp = 0; pp < C->nIPairs; pp++)
  { 
    size_t kk = P[pp*2];
    size_t ll = P[pp*2+1];
    {
      double d = eudist3(X+3*kk, X+3*ll);
      // TODO: min/max to avoid branch?
      if(d > C->dInteraction)
      {
        errInt += pow(d - C->dInteraction, 2);
      }
    }
  }

  // Repulsion
  double errVol = 0;

  errVol = errRepulsion(X, nX, 2*C->r0);

  return C->kInt*errInt + C->kVol*errVol + C->kSph*errSph + C->kRad*errRad;
}

double err2(double * X, size_t nX, double * R, uint32_t * P, conf * C )
{
  /* Alternative version with a list of pairs in contact instead of A */

  // Wanted radii
  double errRad = 0;
  if( (R != NULL) & (C->kRad > 0))
  {
    for(size_t kk = 0; kk<nX; kk++)
    {
      if(isfinite(R[kk]) == 1)
      {
        const double r = norm3(X+kk*3);
        errRad += pow(r-R[kk], 2);
      }
    }
    //  printf("errRad: %f\n", errRad);
  }


  // Keep inside sphere
  double errSph = 0;
  for(size_t kk = 0 ; kk<nX; kk++)
  {
    const double r = norm3(X+kk*3);
    if( r > 1-C->r0)
    {
      errSph += pow(r + C->r0 -1, 2);
    }
  }

  // Wanted contacts/interactions
  double errInt = 0;
  for(size_t pp = 0; pp < C->nIPairs; pp++)
  { 
    size_t kk = P[pp*2];
    size_t ll = P[pp*2+1];
    {
      double d = eudist3(X+3*kk, X+3*ll);
      if(d > C->dInteraction)
      {
        errInt += pow(d - C->dInteraction, 2);
      }
    }
  }

  // Repulsion
  double errVol = 0;
  for(size_t kk = 0; kk < nX; kk++)
  { 
    for(size_t ll = kk+1; ll < nX; ll++)
    {            
      if(fabs(X[3*kk] - X[3*ll]) < 2*C->r0) // This line doubles the speed!
      {
        double d = eudist3(X+3*kk, X+3*ll);
        if( d < 2*C->r0)
        {
          errVol += pow(d - 2*C->r0, 2);
        }
      }
    }
  }

  return C->kInt*errInt + C->kVol*errVol + C->kSph*errSph + C->kRad*errRad;
}


double err(double * X, size_t nX, double * R, uint8_t * A, conf * C )
{

  // Wanted radii
  double errRad = 0;
  if( (R != NULL) & (C->kRad > 0))
  {
    for(size_t kk = 0; kk<nX; kk++)
    {
      if(isfinite(R[kk]) == 1)
      {
        const double r = norm3(X+kk*3);
        errRad += pow(r-R[kk], 2);
      }
    }
    //    printf("errRad: %f\n", errRad);
  }

  // Keep inside sphere
  double errSph = 0;
  for(size_t kk = 0 ; kk<nX; kk++)
  {
    const double r = norm3(X+kk*3);
    if( r > 1-C->r0)
    {
      errSph += pow(r + C->r0 -1, 2);
    }
  }

  // Wanted contacts/interactions
  double errInt = 0;
  for(size_t kk = 0; kk < nX; kk++)
  { 
    for(size_t ll = kk+1; ll < nX; ll++)
    {
      if(A[kk + nX*ll] == 1)
      {
        double d = eudist3(X+3*kk, X+3*ll);
        if(d > C->dInteraction)
        {
          errInt += pow(d - C->dInteraction, 2);
        }
      }
    }
  }

  // Repulsion
  double errVol = 0;
  for(size_t kk = 0; kk < nX; kk++)
  { 
    for(size_t ll = kk+1; ll < nX; ll++)
    {            
      double d = eudist3(X+3*kk, X+3*ll);
      if( d < 2*C->r0)
      {
        errVol += pow(d - 2*C->r0, 2);
      }
    }
  }

  return C->kInt*errInt + C->kVol*errVol + C->kSph*errSph + C->kRad*errRad;
}

void grad(double * X, size_t nX, double * R, uint8_t * A, double * G, conf * C)
{
  // Reset G
  for(size_t kk = 0; kk<nX*3; kk++)
    G[kk] = 0;

  // Radial positioning
  if(C->kRad > 0)
  {
    for(size_t kk = 0; kk<nX; kk++)
    {
      if(isfinite(R[kk]) == 1)
      {
        double r = norm3(X+kk*3);
        double re = 0;
        if(r > 0)
          re = 2*1/r*(r-R[kk]);
        for(int idx = 0; idx<3; idx++)
        {
          G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
        }
      }
    }
  }

  // Keep inside sphere
  for(size_t kk = 0; kk<nX; kk++)
  {
    double r = norm3(X+kk*3);
    if(r > 1-C->r0)
    {
      double re = 2 / r * (r-(1-C->r0));
      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kSph*X[kk*3+idx]*re;
      }
    }
  }

  // Wanted interactions
  for(size_t kk = 0; kk < nX; kk++)
  { 
    for(size_t ll = kk+1; ll < nX; ll++)
    {
      if(A[kk + nX*ll] == 1)
      {
        double d = eudist3(X+3*kk, X+3*ll);
        if(d > C->dInteraction)
        {
          for(int idx = 0; idx<3; idx++)
          {
            G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
            G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
          }
        }
      }
    }
  }

  // Repulsion
  for(size_t kk = 0; kk < nX; kk++)
  { 
    for(size_t ll = kk+1; ll < nX; ll++)
    {
      double d = eudist3(X+3*kk, X+3*ll);
      if( d < 2*C->r0)
      {
        for(int idx = 0; idx<3; idx++)
        {
          G[3*kk+idx] += C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
          G[3*ll+idx] -= C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
        }
      }
    }
  }

  return;
}

void grad2(double * X, size_t nX, double * R, uint32_t * I, double * G, conf * C)
{
  // Reset G
  for(size_t kk = 0; kk<nX*3; kk++)
    G[kk] = 0;

  // Radial positioning
  if(C->kRad > 0)
  {
    for(size_t kk = 0; kk<nX; kk++)
    {
      if(isfinite(G[kk]) == 1)
      {
        double r = norm3(X+kk*3);
        double re = 0;
        if(r > 0)
          re = 2*1/r*(r-R[kk]);
        for(int idx = 0; idx<3; idx++)
        {
          G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
        }
      }
    }
  }



  // Keep inside sphere
  for(size_t kk = 0; kk<nX; kk++)
  {
    double r = norm3(X+kk*3);
    if(r > 1-C->r0)
    {
      double re = 2 / r * (r-(1-C->r0));
      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kSph*X[kk*3+idx]*re;
      }
    }
  }

  // Wanted interactions
  for(size_t pp = 0; pp < C->nIPairs; pp++)
  { 
    size_t kk = I[pp*2];
    size_t ll = I[pp*2+1];

    double d = eudist3(X+3*kk, X+3*ll);
    if(d > C->dInteraction)
    {
      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
        G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
      }
    }
  }

  // Repulsion
  for(size_t kk = 0; kk < nX; kk++)
  { 
    for(size_t ll = kk+1; ll < nX; ll++)
    {
      double d = eudist3(X+3*kk, X+3*ll);
      if( d < 2*C->r0)
      {
        for(int idx = 0; idx<3; idx++)
        {
          G[3*kk+idx] += C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
          G[3*ll+idx] -= C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
        }
      }
    }
  }

  return;
}

void grad3(double * restrict X, 
    const size_t nX, 
    double * restrict R, 
    uint32_t * restrict I, 
    double * restrict G, 
    const conf * restrict C)
{
  // Reset G
  for(size_t kk = 0; kk<nX*3; kk++)
    G[kk] = 0;

  // Radial positioning
  if(C->kRad > 0)
  {
    for(size_t kk = 0; kk<nX; kk++)
    {
      if(isfinite(R[kk]) == 1)
      {
        double r = norm3(X+kk*3);
        double re = 0;
        // if(r > 0)
        re = 2*1/r*(r-R[kk]);
        for(int idx = 0; idx<3; idx++)
        {
          G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
        }
      }
    }
  }

  // Keep inside sphere
  for(size_t kk = 0; kk<nX; kk++)
  {
    double r = norm3(X+kk*3);
    if(r > 1-C->r0)
    {
      double re = 2 / r * (r-(1-C->r0));
      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kSph*X[kk*3+idx]*re;
      }
    }
  }

  // Wanted interactions
  for(size_t pp = 0; pp < C->nIPairs; pp++)
  { 
    size_t kk = I[pp*2];
    size_t ll = I[pp*2+1];

    double d = eudist3(X+3*kk, X+3*ll);
    assert(kk != ll);
#ifndef NDEBUG
if( !(d>0) )
{
  printf("Strange distance between interacting points!\n");
  printf("@%p : %f %f %f\n", X+3*kk, X[3*kk], X[3*kk+1], X[3*kk+2]);
  printf("@%p : %f %f %f\n", X+3*ll, X[3*ll], X[3*ll+1], X[3*ll+2]);
  exit(1);
}
#endif
    if(d > C->dInteraction)
    {
      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
        G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
      }
    }
  }

  // Repulsion
  gradRepulsion(X, G, nX, 2*C->r0, C->kVol);

  return;
}

void grad4(double * restrict X, 
    const size_t nX, 
    double * restrict R, 
    uint32_t * restrict I, 
    double * restrict G, 
    const conf * restrict C)
{
  // Reset G
  for(size_t kk = 0; kk<nX*3; kk++)
    G[kk] = 0;

  // Radial positioning
  if(C->kRad > 0)
  {
    for(size_t kk = 0; kk<nX; kk++)
    {
      if(isfinite(R[kk]) == 1)
      {


      double r = norm3(X+kk*3);
      double re = 0;
      // if(r > 0)
      re = 2*1/r*(r-R[kk]);
      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
      }
      }
    }
  }

  // Keep inside sphere
  const double rmax = 1-C->r0;
  for(size_t kk = 0; kk<nX; kk++)
  {
    double r = norm3(X+kk*3);
    if(r > rmax)
    {
      double re = 2 / r * (r-(1-C->r0));
      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kSph*X[kk*3+idx]*re;
      }
    }
  }

  // Wanted interactions
  const double dInteraction = C->dInteraction;
  const double dInteraction2 = pow(dInteraction, 2);

  for(size_t pp = 0; pp < C->nIPairs; pp++)
  { 
    size_t kk = I[pp*2];
    size_t ll = I[pp*2+1];

    double d2 = eudist3sq(X+3*kk, X+3*ll);
    assert(kk != ll);
    assert(d2>0);
    if(d2 > dInteraction2)
    {
      double d= sqrt(d2);

      for(int idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - dInteraction);
        G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - dInteraction);
      }
    }
  }

  // Repulsion
  gradRepulsion(X, G, nX, 2*C->r0, C->kVol);

  return;
}
