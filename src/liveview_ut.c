#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include "liveview.h"
#include <pthread.h>
#include <unistd.h>

int main(int argc, char ** argv)
{
    if(argc > 1)
    {
        printf("WARNING: %s takes no parameters\n", argv[0]);
    }

  size_t N = 3030;

  double *X = malloc(3*N*sizeof(double));
  assert(X != NULL);
  uint8_t * L = malloc(N*sizeof(uint8_t));
  assert(L != NULL);


  for(size_t kk = 0; kk<3*N ; kk++)
  {
    X[kk] = 2*(.5-((double) rand() / (double) RAND_MAX));
  }
  for(size_t kk = 0; kk<N ; kk++)
  {
    L[kk] = kk % 23;
  }


  liveXLview xlview;
  xlview.X = X;
  xlview.L = L;
  xlview.N = N;
  xlview.quit = 0;

  pthread_t th;
  pthread_create(&th, // thread
      NULL, // pthread_attrib_t
      liveview_t, // function
      &xlview); // arg

  for(size_t iter = 0; iter<5000; iter++)
  {
  for(size_t kk = 0; kk<3*N; kk++)
  {
    X[kk] += 0.01*(((double) rand() / (double) RAND_MAX)-.5);
  }
  usleep(100);
  }

  xlview.quit = 1;
  pthread_join(th, NULL);
  // liveview(X, L, N);



  free(X);
  free(L);
  return 0;
}
