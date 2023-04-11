#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "hsvrgb.c"
#include "hsvrgb.h"

static void usage(char * cmd)
{
  fprintf(stdout, "%s r g b\n", cmd);
  return;
}

double static eudist3(const double * A, const double * B)
{
  /* Euclidean distance between two 3D-vectors */
  return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}

int main(int argc, char ** argv)
{


  double * rgb = malloc(3*sizeof(double));
  double * rgb2 = malloc(3*sizeof(double));
  double * hsv = malloc(3*sizeof(double));
  double * hsv2 = malloc(3*sizeof(double));
  double * rgb3 = malloc(3*sizeof(double));

  if(argc == 4)
  {
    for(int idx = 0; idx<3; idx++)
    {
      rgb[idx] = hsv[idx] = atof(argv[idx+1]);
    }

    rgb2hsv(rgb, hsv2);
    hsv2rgb(hsv, rgb2);

    fprintf(stdout, "as RGB: %f %f %f\n", rgb[0], rgb[1], rgb[2]);
    fprintf(stdout, " -> HSV: %f %f %f\n", hsv2[0], hsv2[1], hsv2[2]);

    fprintf(stdout, "as HSV: %f %f %f\n", hsv[0], hsv[1], hsv[2]);
    fprintf(stdout, " -> RGB: %f %f %f\n", rgb2[0], rgb2[1], rgb2[2]);

    hsv2rgb(hsv2, rgb3);

    fprintf(stdout, "Round trip, RGB->HSV->RGB\n");
    fprintf(stdout, "% f, % f, % f\n", 
        rgb[0] - rgb3[0],
        rgb[1] - rgb3[1],
        rgb[2] - rgb3[2]);
    if(eudist3(rgb, rgb3)>0.01)
    {
      fprintf(stderr, "That's wrong!\n");
      exit(-1);
    }
  }

  double delta = 0.001;
  fprintf(stdout, "Testing rgb->hsv->rgb with step %f\n", delta);
  /* Can't test the other way around since hsv values are not unique,
   * i.e., (0,1,1) and (1,1,1) are the same 
   */
  size_t nTests = 0;
  for(double r = 0; r<=1; r+=delta)
  {
    for(double g = 0; g<=1; g+=delta)
    {
      for(double b = 0; b<=1; b+=delta)
      {
        rgb[0] = r; rgb[1] = g; rgb[2] = b;
        rgb2hsv(rgb, hsv);
        assert(rgb[0] == r);
        assert(rgb[1] == g);
        assert(rgb[2] == b);

        double h = hsv[0]; double s = hsv[1]; double v = hsv[2];
        hsv2rgb(hsv, rgb2);
        assert(hsv[0] == h);
        assert(hsv[1] == s);
        assert(hsv[2] == v);

        if(eudist3(rgb, rgb2) > 0.01)
        {
          printf("Error!\n");
          printf("rgb: %f %f %f\n", rgb[0], rgb[1], rgb[2]);
          printf("hsv: %f %f %f\n", hsv[0], hsv[1], hsv[2]);
          printf("rgb2: %f %f %f\n", rgb2[0], rgb2[1], rgb2[2]);
          exit(-1);
        }
        nTests++;
      }
    }
  }
  fprintf(stdout, "%zu tests passed!\n", nTests);

  free(rgb);
  free(rgb2);
  free(rgb3);
  free(hsv);
  free(hsv2);

  return 0;
}
