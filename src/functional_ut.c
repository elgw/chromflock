#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>

#include "functional.h"

static double clockdiff(struct timespec* start, struct timespec * finish)
{
  double elapsed = (finish->tv_sec - start->tv_sec);
  elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
  return elapsed;
}


void test_ellipsoid(size_t N, // number of points
    double * X, // 3xN coordinates of points
    size_t nPairs, // number of pair constraints
    uint32_t * P, // Their 2xnPairs indices
    double * R, // wanted radius
    double * G1, double * G2) //storage for gradeints
{
  printf("--> Testing _ellipsoidal_ geometry\n");
  printf("  --> Spherical ellipsoid vs sphere\n");

  conf fconf;
  fconf.r0 = 0.03;
  fconf.kVol = 0; // volume exclusion -- repulsion
  fconf.kDom = 0;
  fconf.kInt = 0;
  fconf.kRad = 0;
  fconf.dInteraction = 2.1*fconf.r0;
  fconf.nIPairs = nPairs;
  fconf.E = NULL;

  for(int kVol = 0; kVol<2; kVol++)
  {
  for(int kInt = 0; kInt<2; kInt++)
  {
    for(int kRad = 0; kRad<2; kRad++)
    {
      for(int kDom = 0; kDom<2; kDom++)
      {
        fconf.r0 = 0.03;
        fconf.kVol = kVol; // volume exclusion -- repulsion
        fconf.kDom = kDom;
        fconf.kInt = kInt;
        fconf.kRad = kRad;
        fconf.dInteraction = 2.1*fconf.r0;
        fconf.nIPairs = nPairs;
        fconf.E = NULL;

        //        printf("--> Spherical\n");
        double err = err3(X, N, R, P, &fconf);
        //        printf("err: %.10f\n", err);
        grad3(X, N, R, P, G1, &fconf);

        //      printf("--> Spherical ellipsoid\n");
        double r = 1;

        fconf.E = elli_new(r, r, r);
        fconf.Es = elli_new(r-fconf.r0, r-fconf.r0, r-fconf.r0);

        double erre = err3(X, N, R, P, &fconf);
        //      printf("erre: %.10f\n", erre);
        grad3(X, N, R, P, G2, &fconf);

        if( fabs(err-erre) > 1e-9)
        {
          printf("kVol: %d, kInt: %d, kRad: %d, kDom: %d\n", kVol, kInt, kRad, kDom);
          printf("err: %f, erre: %f\n", err, erre);
          assert(0);
        }

        //    printf("--> Verifying equal gradient\n");
        for(size_t kk = 0; kk<N; kk++)
        {
        for(size_t ll = 0; ll<3; ll++)
        {
          if(fabs(G1[3*kk+ll]-G2[3*kk+ll])>1e-5)
          {
            printf("X = [%f, %f, %f]\n", X[3*kk], X[3*kk+1], X[3*kk+2]);
            printf("Gs(%zu, %zu)=Gs[%zu] = %.10f\n", kk, ll, 3*kk+ll, G1[kk]);
            printf("Ge(%zu, %zu)=Ge[%zu] = %.10f\n", kk, ll, 3*kk+ll, G2[kk]);
            assert(0);
          }
        }
        }
      }
    }
  }
  }

  double a = 1; double b = 0.8; double c = 0.8;
  fconf.r0 = 0.03;
  fconf.E = elli_new(a, b, c);
  fconf.Es = elli_new(a-fconf.r0, b-fconf.r0, c-fconf.r0);
  fconf.dInteraction = 2.1*fconf.r0;
  fconf.nIPairs = nPairs;

  double delta = 1e-9;

  printf("  --> Numerical gradient\n");
  for(int kVol = 0; kVol<2; kVol++) {
    for(int kDom = 0; kDom<2; kDom++) {
      for(int kInt = 0; kInt<2; kInt++) {
        for(int kRad = 0; kRad<2; kRad++) {

          fconf.kVol = kVol; // volume exclusion -- repulsion
          fconf.kDom = kDom;
          fconf.kInt = kInt;
          fconf.kRad = kRad;

          grad3(X, N, R, P, G1, &fconf);

          double e0 = err3(X, N, R, P, &fconf);
          int gotError = 0;
          for(size_t kk = 0; kk<3*N; kk++)
          {
            double saved = X[kk];
            X[kk] = X[kk] + delta;
            double e1 = err3(X, N, R, P, &fconf);
            double numd = (e1-e0)/delta;
            double diff = fabs(G1[kk]-numd);
            if(diff > 1e-3)
            {
              printf("%f, (%f) ", R[kk], saved);
            printf("kVol: %d kDom: %d, kInt: %d, kRad: %d\n", kVol, kDom, kInt, kRad);
              printf("de/dx_%zu = %f, num = %f, delta= %f, q = %f) \n",
                  kk, G1[kk], numd, diff, G1[kk]/numd);
              gotError = 1;
              getchar();

            }
            X[kk] = saved;
          }

          if(gotError == 1)
          {
            printf("Numerical and analytical gradient does not match!\n");
            printf("kVol: %d kDom: %d, kInt: %d, kRad: %d\n", kVol, kDom, kInt, kRad);
            exit(-1);
          }
        }
      }
    }
  }

  printf("/// Testing ellipsoidal geometry\n");

}

int main(int argc, char ** argv)
{
  /* Unit test.
   * Compare gradient to numerical gradient */


  size_t N = 600;
  int rep = 10000;

  if(argc > 1)
  {
    N = atol(argv[1]);
  }
  if(argc > 2)
  {
      rep = atoi(argv[2]);
  }

  printf("Using %zu beads\n", N);

  double * X = malloc(3*N*sizeof(double));
  assert(X != NULL);
  double * R = malloc(N*sizeof(double));
  assert(R != NULL);
  double * G1 = malloc(3*N*sizeof(double));
  assert(G1 != NULL);
  double * G2 = malloc(3*N*sizeof(double));
  assert(G2 != NULL);
  double * G3 = malloc(3*N*sizeof(double));
  assert(G3 != NULL);

  printf("Constructing contact matrix (TODO: use only pair list)\n");
  uint8_t * A = malloc(N*N*sizeof(uint8_t));
  assert(A != NULL);

  for(size_t kk = 0; kk<N*N; kk++)
  {
    A[kk] = 0;
  }
  for(size_t kk = 1; kk<N*N; kk=kk+N+1)
  {
    A[kk] = 1;
    A[kk+N-1] = 1;
  }

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
  assert(P!=NULL);
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
    X[kk] = 2.1*(0.5-rand()/(double) RAND_MAX);

  double delta = 1e-8;
  struct timespec ta, tb;
  conf fconf;

  /* TESTS */



  printf("--> Timing the fastest versions using %d repetitions\n", rep);
  {
    fconf.r0 = 0.03;
    fconf.kVol = 1; // volume exclusion -- repulsion
    fconf.kDom = 1;
    fconf.kInt = 1;
    fconf.kRad = 1;
    fconf.dInteraction = 2.1*fconf.r0;
    fconf.nIPairs = nPairs;
    fconf.E = NULL;



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
    for(int kDom = 0; kDom<2; kDom++) {
      for(int kInt = 0; kInt<2; kInt++) {
        for(int kRad = 0; kRad<2; kRad++) {

          fconf.r0 = 0.03;
          fconf.kVol = kVol; // volume exclusion -- repulsion
          fconf.kDom = kDom;
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
            printf("kVol: %d kDom: %d, kInt: %d, kRad: %d\n", kVol, kDom, kInt, kRad);
            exit(-1);
          }
        }
      }
    }
  }

  printf("\n");

  test_ellipsoid(N, X, nPairs, P,  R, G1, G2);


  fprintf(stdout, "Testing gradients\n");
  fprintf(stdout, "Testing all combinations of forces\n");

  size_t nerrors = 0;

  for(int kVol = 0; kVol<2; kVol++) {
    for(int kDom = 0; kDom<2; kDom++) {
      for(int kInt = 0; kInt<2; kInt++) {
        for(int kRad = 0; kRad<2; kRad++) {

          fconf.r0 = 0.03;
          fconf.kVol = kVol; // volume exclusion -- repulsion
          fconf.kDom = kDom;
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
    for(int kDom = 0; kDom<2; kDom++) {
      for(int kInt = 0; kInt<2; kInt++) {
        for(int kRad = 0; kRad<2; kRad++) {

          fconf.r0 = 0.03;
          fconf.kVol = kVol; // volume exclusion -- repulsion
          fconf.kDom = kDom;
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
            printf("kVol: %d kDom: %d, kInt: %d, kRad: %d\n", kVol, kDom, kInt, kRad);
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
