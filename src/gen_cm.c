#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

void usage(char * binName)
{
  fprintf(stdout, "Generates a contact probability matrix, A, for testing\n");
  fprintf(stdout, "Usage:\n");
  fprintf(stdout, "%s file.dat N M\n", binName);
  fprintf(stdout, "Where N is the number of beads and M the number of contacts\n");
  return;
}

int main(int argc, char ** argv)
{

  if(argc<2)
  {
    usage(argv[0]);
    exit(-1);
  }


  FILE * fout = fopen(argv[1], "w");
  
  if(fout == NULL)
  {
    fprintf(stderr, "Could not open %s for writing\n", argv[1]);
    exit(-1);
  }

  int offDiagonal = 1;
  size_t N = 100;
  size_t nrand = 1; // number of random contacts

  if(argc>2)
  {
    N = atol(argv[2]);
  }

  if(argc>3)
  {
    offDiagonal = atoi(argv[3]);
  }
  fprintf(stdout, "Settings:\n");
  fprintf(stdout, "N: %zu\n", N);
  fprintf(stdout, "1-diagonal: %d\n", offDiagonal);
  fprintf(stdout, "N rand: %zu\n", nrand);

  double * A = malloc(N*N*sizeof(double));

  memset(A, 0, N*N*sizeof(double));


  if(offDiagonal)
  {
    
    for(size_t kk = 1 ; kk < N*N ; kk=kk+N+1)
    {
      A[kk] = 1;
    }
 for(size_t kk = N ; kk < N*N ; kk=kk+N+1)
    {
      A[kk] = 1;
    }
  }

  // Random contacts
  size_t randset = 0;
  while(randset < nrand)
  {
    size_t kk = N*rand()/RAND_MAX;
    size_t ll = N*rand()/RAND_MAX;
    assert(kk<N);
    assert(ll<N);
    if(A[kk + N*ll] == 0)
    {
      randset++;
      A[kk + N*ll] = 1;
      A[ll + N*kk] = 1;
    }
  }


  // Count zeros and ones
  size_t ones = 0;
  size_t zeros = 0;
  for(size_t kk = 0; kk<N*N; kk++)
  {
    if(A[kk] == 1)
      ones++;
    if(A[kk] == 0)
      zeros++;
  }

  // Check symmetry
  for(size_t kk = 0; kk < N; kk++)
  {
    for(size_t ll = 0; ll < N; ll++)
    {
      assert( A[kk + N*ll] == A[ll + N*kk]);
    }
  }

  fprintf(stdout, "%zu zeros, %zu ones\n", zeros, ones);

  size_t nwritten = fwrite(A, sizeof(double), N*N, fout);
  fprintf(stdout, "Wrote %zu bytes / %zu elements\n", nwritten*sizeof(double), N*N);

  fclose(fout);
}


