#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>

#include "functional.h"

double static clockdiff(struct timespec* start, struct timespec * finish)
{
  double elapsed = (finish->tv_sec - start->tv_sec);
  elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
  return elapsed;
}

int main(int argc, char ** argv)
{
  /* Unit test. 
   * Compare gradient to numerical gradient */


  size_t N = 600;
  if(argc > 1)
  {
    N = atol(argv[1]);
  }

  if(argc>1)
  {
    N = atol(argv[1]);
  }

  double * X = malloc(3*N*sizeof(double));
  double * R = malloc(N*sizeof(double));
  double * G1 = malloc(3*N*sizeof(double));
  double * G2 = malloc(3*N*sizeof(double));
  double * G3 = malloc(3*N*sizeof(double));
  uint8_t * A = malloc(N*N*sizeof(uint8_t));

  for(size_t kk = 0; kk<N*N; kk++)
  {
    A[kk] = 0;
  }
  for(size_t kk = 1; kk<N*N; kk=kk+N+1)
  {
    A[kk] = 1;
    A[kk+N-1] = 1;
  }

  /* Construct list of pairs */
  // Count them
  size_t nPairs = 0;
  for(size_t kk = 0; kk<N; kk++)
    for(size_t ll = kk+1; ll<N; ll++)
    {
      if(A[kk+ll*N] == 1)
        nPairs++;
    }

  printf("Found %zu pairs\n", nPairs);
  // Construct the list
  size_t writepos = 0;
  uint32_t * P = malloc(nPairs*2*sizeof(uint32_t));
  for(size_t kk = 0; kk<N; kk++)
    for(size_t ll = kk+1; ll<N; ll++)
    {
      if(A[kk+ll*N] == 1)
      {
        P[writepos++] = kk;
        P[writepos++] = ll;
      }
    }
  assert(writepos == 2*nPairs);



  for(size_t kk = 0; kk<N; kk++)
    R[kk] = rand()/(double) RAND_MAX;

  for(size_t kk = 0; kk<N*3; kk++)
    X[kk] = 2*(0.5-rand()/(double) RAND_MAX);

  double delta = 1e-8;
  struct timespec ta, tb;
  conf fconf;

  /* TESTS */


  /* Timings just for the fastest versions */
  {
  fconf.r0 = 0.03;
  fconf.kVol = 1; // volume exclusion -- repulsion
  fconf.kSph = 1;
  fconf.kInt = 1;
  fconf.kRad = 1;
  fconf.dInteraction = 2.1*fconf.r0;
  fconf.nIPairs = nPairs;

  const int rep = 10000;

  clock_gettime(CLOCK_MONOTONIC, &ta);
  for(int kk = 0; kk<rep; kk++)
  {
    grad3(X, N, R, P, G3, &fconf);
  }
  clock_gettime(CLOCK_MONOTONIC, &tb);
  printf("   Analytic grad3 took: %f s\n", clockdiff(&ta, &tb)/rep);

  clock_gettime(CLOCK_MONOTONIC, &ta);
  for(int kk = 0; kk<rep; kk++)
  {
    grad4(X, N, R, P, G3, &fconf);
  }
  clock_gettime(CLOCK_MONOTONIC, &tb);
  printf("   Analytic grad4 took: %f s\n", clockdiff(&ta, &tb)/rep);

  clock_gettime(CLOCK_MONOTONIC, &ta);
  double e3 = 0;
  for(int kk = 0; kk<rep; kk++)
  {
  e3 += err3(X, N, R, P, &fconf);
  }
  clock_gettime(CLOCK_MONOTONIC, &tb);
  printf("   Analytic err3 took: %f s\n", clockdiff(&ta, &tb)/rep);
  }

  /* Test that it works. Note N=1000 takes about 1 min */
  printf("Comparing analytic with numerical gradient\n");
  printf("Using combinations of the parts of the functionals\n");

  for(int kVol = 0; kVol<2; kVol++) {
    for(int kSph = 0; kSph<2; kSph++) {
      for(int kInt = 0; kInt<2; kInt++) {
        for(int kRad = 0; kRad<2; kRad++) {

          fconf.r0 = 0.03;
          fconf.kVol = kVol; // volume exclusion -- repulsion
          fconf.kSph = kSph;
          fconf.kInt = kInt;
          fconf.kRad = kRad;
          fconf.dInteraction = 2.1*fconf.r0;
          fconf.nIPairs = nPairs;

          grad3(X, N, R, P, G3, &fconf);

          double e0 = err(X, N, R, A, &fconf);
          int gotError = 0;
          for(size_t kk = 0; kk<3*N; kk++)
          {
            double saved = X[kk];
            X[kk] = X[kk] + delta;
            double e1 = err(X, N, R, A, &fconf);
            double numd = (e1-e0)/delta;
            double diff = fabs(G3[kk]-numd);
            if(diff > 1e-4)
            {
              printf("de/dx_%zu = %f, num = %f, delta= %f, q = %f) \n", 
                  kk, G3[kk], numd, diff, G3[kk]/numd);
              gotError = 1;

            }
            X[kk] = saved;
          }

          if(gotError == 1)
          {
            printf("Numerical and analytical gradient does not match!\n");
            printf("kVol: %d kSph: %d, kInt: %d, kRad: %d\n", kVol, kSph, kInt, kRad);
            exit(-1);
          }


        }}}}

  printf("\n");


  fprintf(stdout, "Testing gradients\n");
  fprintf(stdout, "Testing all combinations of forces\n");

  size_t nerrors = 0;

  for(int kVol = 0; kVol<2; kVol++) {
    for(int kSph = 0; kSph<2; kSph++) {
      for(int kInt = 0; kInt<2; kInt++) {
        for(int kRad = 0; kRad<2; kRad++) {

          fconf.r0 = 0.03;
          fconf.kVol = kVol; // volume exclusion -- repulsion
          fconf.kSph = kSph;
          fconf.kInt = kInt;
          fconf.kRad = kRad;
          fconf.dInteraction = 2*fconf.r0;
          fconf.nIPairs = nPairs;


          clock_gettime(CLOCK_MONOTONIC, &ta);
          grad3(X, N, R, P, G3, &fconf);
          clock_gettime(CLOCK_MONOTONIC, &tb);
          printf("   Analytic grad3 took: %f s\n", clockdiff(&ta, &tb));

          clock_gettime(CLOCK_MONOTONIC, &ta);
          grad(X, N, R, A, G1, &fconf);
          clock_gettime(CLOCK_MONOTONIC, &tb);
          printf("   Analytic grad  took: %f s\n", clockdiff(&ta, &tb));

          clock_gettime(CLOCK_MONOTONIC, &ta);
          grad(X, N, R, A, G2, &fconf);
          clock_gettime(CLOCK_MONOTONIC, &tb);
          printf("   Analytic grad2 took: %f s\n", clockdiff(&ta, &tb));

          printf("  G1 vs G2\n");

          for(size_t kk = 0; kk < N*3; kk++)
          {
            if(fabs(G1[kk]-G2[kk]) > 1e-7)
              //    if(fabs(G1[kk]) > 1e-7)
            {
              printf("G1[%zu]=%f G3[%zu]=%f\n", kk, G1[kk], kk, G3[kk]);
              nerrors++;
            }
          }
          printf("  G1 vs G3\n");

          for(size_t kk = 0; kk < N*3; kk++)
          {
            if(fabs(G1[kk]-G3[kk]) > 1e-7)
              //    if(fabs(G1[kk]) > 1e-7)
            {
              printf("G1[%zu]=%f G3[%zu]=%f\n", kk, G1[kk], kk, G3[kk]);
              nerrors++;
            }
          }

          if(nerrors > 0)
          {
            printf("There were errors !\n");
            exit(-1);
          }

        }}}}


  printf("\n");
  printf("Comparing error functions\n");

  for(int kVol = 0; kVol<2; kVol++) {
    for(int kSph = 0; kSph<2; kSph++) {
      for(int kInt = 0; kInt<2; kInt++) {
        for(int kRad = 0; kRad<2; kRad++) {

          fconf.r0 = 0.03;
          fconf.kVol = kVol; // volume exclusion -- repulsion
          fconf.kSph = kSph;
          fconf.kInt = kInt;
          fconf.kRad = kRad;
          fconf.dInteraction = 2*fconf.r0;
          fconf.nIPairs = nPairs;

          clock_gettime(CLOCK_MONOTONIC, &ta);
          double e0 = err(X, N, R, A, &fconf);
          clock_gettime(CLOCK_MONOTONIC, &tb);
          printf("   E0 took: %f s\n", clockdiff(&ta, &tb));

          clock_gettime(CLOCK_MONOTONIC, &ta);
          double e2 = err2(X, N, R, P, &fconf);
          clock_gettime(CLOCK_MONOTONIC, &tb);
          printf("   E2 took: %f s\n", clockdiff(&ta, &tb));

          clock_gettime(CLOCK_MONOTONIC, &ta);
          double e3 = err3(X, N, R, P, &fconf);
          clock_gettime(CLOCK_MONOTONIC, &tb);
          printf("   E3 took: %f s\n", clockdiff(&ta, &tb));

          if(fabs(e0-e2)>1e-6 || fabs(e0-e3) > 1e-6)
          {
            printf("Error functions give different results!\n");
            printf("kVol: %d kSph: %d, kInt: %d, kRad: %d\n", kVol, kSph, kInt, kRad);
            printf("err: %f err2: %f err3: %f\n", e0, e2, e3);
            exit(1);
          }

        }}}}

  printf("\n");
  printf("All tests passed!\n");

  free(G1);
  free(G2);
  free(G3);
  free(R);
  free(X);
  free(A);
  free(P);

  return 0;
}
