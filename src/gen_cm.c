#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

void usage(char * binName)
{
  fprintf(stdout, "Generates a contact probability matrix, A, for testing\n");
  fprintf(stdout, "Usage:\n");
  fprintf(stdout, "%s name N M\n", binName);
  fprintf(stdout, " name: output name, A_name.double and L_name.uint8 will be created.\n");
  fprintf(stdout, " N:  the number of beads.\n");
  fprintf(stdout, " M:  the number of contacts\n");

  return;
}

int main(int argc, char ** argv)
{

  if(argc<2)
  {
    usage(argv[0]);
    exit(-1);
  }

  char * name = argv[1];
  char * aFile = malloc(1024*sizeof(char));
  char * lFile = malloc(1024*sizeof(char));

  sprintf(aFile, "A_%s.double", name);
  sprintf(lFile, "L_%s.uint8", name);

  FILE * aOut = fopen(aFile, "w");
  FILE * lOut = fopen(lFile, "w");
  
  if(aOut == NULL)
  {
    fprintf(stderr, "Could not open %s for writing\n", aFile);
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

  size_t nwritten = fwrite(A, sizeof(double), N*N, aOut);
  fprintf(stdout, "Wrote %zu bytes / %zu elements to %s\n", nwritten*sizeof(double), N*N, aFile);
  fclose(aOut);

  uint8_t * L = malloc(N*sizeof(uint8_t));
  for(size_t kk = 0; kk<N; kk++)
  {
    L[kk] = 1;
  }
  
  nwritten = fwrite(L, sizeof(uint8_t), N, lOut);
  fprintf(stdout, "Wrote %zu bytes / %zu elements to %s\n", nwritten*sizeof(uint8_t), N, lFile);
  fclose(lOut);
}


