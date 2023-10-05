#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "functional.h"
#define INLINED inline __attribute__((always_inline))

INLINED static double norm3(double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
}

int main(int argc, char ** argv)
{
    if(argc > 1)
    {
        printf("WARNING: %s takes no arguments\n", argv[0]);
    }
    size_t N = 33000;
    double vq = 0.2;
    double Vd = 4/3*M_PI;
    double d = cbrt( 3.0*vq*Vd / (4.0*N*M_PI) );

    printf("Testing with %zu points and d=%f\n", N, d);
    double * X = malloc(3*N*sizeof(double));
    assert(X != NULL);

    // 1. Create N beads in a sphere
    for(size_t kk = 0; kk< N; kk++)
    {
        int accepted = 0;
        while(accepted == 0)
        {
            for(int idx =0; idx<3; idx++)
            {
                X[3*kk+idx] = 2*(rand()/(double) RAND_MAX-.5);
            }
            if(norm3(X+3*kk) <= 1)
            { accepted = 1;}
        }
    }

    size_t nIter = 20000000/N;
    printf("Using %zu iterations\n", nIter);
    double error = 0;
    for(size_t kk = 0; kk<nIter; kk++)
    {
        error += errRepulsion(X, N , d);
    }
    printf("ignore: %f\n", error); // not relevant, just to know that nothing is optimized away

}
