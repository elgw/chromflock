#include "balance.h"

static void rsums0(double * R, double * A, size_t N)
{
  // Sum of all rows
  for(size_t kk = 0; kk<N; kk++)
  {
    R[kk] = 0;
    for(size_t ll = 0; ll<N; ll++)
    {
      double v = A[kk*N + ll];
      R[kk] += v;
    }
  }
}

double balance(double * A, size_t N)
{
  const int verbose = 0;

  double * R = malloc(N*sizeof(double));
  assert(R != NULL);

  for(size_t iter = 0; iter<24; iter++)
  {
    if(verbose) {
    printf("."); fflush(stdout);
    }
    rsums0(R, A, N);
    for(size_t kk = 0; kk<N; kk++)
      R[kk] = sqrt(R[kk]);

    for(size_t kk = 0; kk<N; kk++)
    {
      for(size_t ll = 0; ll<N; ll++)
      {
        double s = R[kk]*R[ll];
        if(s == 0)
        {
          s = 1;
        }
        A[kk+ll*N] /= s;
      }
    }
  }
  if(verbose) {
printf("\n");
  }
  rsums0(R, A, N);
  double maxError = 0;

  /*
   * Check error and ignore empty rows/columns
   */

  size_t nzeros = 0;
  for(size_t kk = 0; kk <  N; kk++)
  {
    if(R[kk] == 0)
    {
      nzeros++;
    }
    else
    {
    double aerror = fabs(R[kk]-1.0);
    if( aerror > maxError)
      maxError = aerror;
    }
  }

  if(nzeros == N)
  {
    printf("   Can't estimate error, all zeros\n");
    maxError = -1;
  } else {
    printf("    Max abs error: %e\n", maxError);
  }

  free(R);

  return maxError;
}
