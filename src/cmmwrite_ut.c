#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <stdint.h>
#include <cairo-svg.h>

#include "cmmwrite.h"
#include "cmmwrite.c"

double static clockdiff(struct timespec* start, struct timespec * finish)
{
  double elapsed = (finish->tv_sec - start->tv_sec);
  elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
  return elapsed;
}

void savecmap()
{
  cairo_t *cr;
  cairo_surface_t *surface;

  surface =
    (cairo_surface_t *)cairo_svg_surface_create("cmm_cmap.svg", 525.0, 525.0);

  cr = cairo_create(surface);


  /* Draw the squares in the background */


  for (int kk=0; kk<5; kk++)
  {
    for (int ll=0; ll<5; ll++)
    {
      int L = ll*5+kk + 1;
      if(L<25)
      {
      double r = (double) cmap[3*L] / 255.0;
      double g = (double) cmap[3*L+1] / 255.0;
      double b = (double) cmap[3*L+2] / 255.0;

      cairo_set_source_rgb(cr, r, g, b);

      int y = 25+kk*100;
      int x = 25+ll*100;

      cairo_rectangle(cr, y, x, 75, 75);
      cairo_fill (cr);
      }
    }
  }

  /* Writing in the foreground */

  char * label = malloc(1024*sizeof(char));

  for (int kk=0; kk<5; kk++)
  {
    for (int ll=0; ll<5; ll++)
    {
      int y = 25+kk*100;
      int x = 25+ll*100-7;

      int L = ll*5+kk + 1;
      if(L<25)
      {

      cairo_set_font_size (cr, 15);

      cairo_select_font_face (cr, "Georgia",
          CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

      cairo_set_source_rgb (cr, 0, 0, 0);

      cairo_move_to(cr, y, x);

      sprintf(label, "chr%d", L);
      if(L==23)
      {
        sprintf(label, "chrX");
      }

      if(L==24)
      {
        sprintf(label, "chrY");
      }

      cairo_show_text(cr, label);
      }
    }
  }


  cairo_destroy (cr);

  cairo_surface_destroy (surface);

  return;
}

int main(int argc, char ** argv)
{

  savecmap();

  size_t nD = 3000;
  size_t nP = 0;

  double * D = malloc(3*nD*sizeof(double));
  uint32_t * P = malloc(2*nP*sizeof(uint32_t));

  for(size_t kk = 0; kk<3*nD; kk++)
  {
    D[kk] = kk;
  }

  uint8_t * L = malloc(nD*sizeof(uint8_t));

  for(size_t kk = 0; kk<nD; kk++)
  {
    L[kk] = (uint8_t) kk%255; // cyclic...
  }


  struct timespec ta, tb;
  clock_gettime(CLOCK_MONOTONIC, &ta);
  cmmwrite("test.cmm", D, nD, 0.1, P, nP, L);
  clock_gettime(CLOCK_MONOTONIC, &tb);

  printf("cmmwrite took %f s\n", clockdiff(&ta, &tb));

  clock_gettime(CLOCK_MONOTONIC, &ta);
  cmmwritez("testz.cmm.gz", D, nD, 0.1, P, nP, L);
  clock_gettime(CLOCK_MONOTONIC, &tb);
  printf("cmmwritez took %f s\n", clockdiff(&ta, &tb));


  free(L);
  free(D);
  free(P);

  return 0;
}
