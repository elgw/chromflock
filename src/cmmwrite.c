/* Example from
https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/volumepathtracer/volumepathtracer.html#markerfiles

<marker_set name="marker set 1">
<marker id="1" x="-6.1267" y="17.44" z="-3.1338"  radius="0.35217"/>
<marker id="2" x="1.5395" y="16.277" z="-3.0339" r="0" g="1" b="1"
radius="0.5" note="An example note"/>
<link id1="2" id2="1" r="1" g="1" b="0" radius="0.17609"/>
</marker_set>
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "cmmwrite.h"

uint8_t cmap[] = {255,255,255, 	
  240,163,255,
  0,117,220, 	
  153,63,0, 	
  76,0,92, 	
  25,25,25, 	
  0,92,49, 	
  43,206,72, 	
  255,204,153, 	
  128,128,128, 	
  148,255,181, 	
  143,124,0, 	
  157,204,0, 	
  194,0,136, 	
  0,51,128, 	
  255,164,5, 	
  255,168,187, 	
  66,102,0, 	
  255,0,16, 	
  94,241,242, 	
  0,153,143, 	
  224,255,102, 	
  116,10,255, 	
  153,0,0, 	
  255,255,128, 	
  255,255,0, 	
  255,80,5};


int cmmwritez(char * fname, double * D, size_t nD, double radius, uint32_t * P, size_t NP, uint8_t * L)
{

  gzFile zf = gzopen(fname, "wb");
  char * line = malloc(1024*sizeof(char));

  if(zf == NULL)
  {
    fprintf(stderr, "Unable to gzopen %s\n", fname);
    return -1;
  }

  sprintf(line, "<marker_set name=\"dump\">\n");
  gzwrite(zf, line, strlen(line));


  for(size_t kk = 0; kk<nD; kk++)
  {

    double r = 1;
    double g = 0;
    double b = 0;

    if(L != NULL)
    {
        uint8_t chr = L[kk];
        chr = chr % 32; // use 33, 34, ... for the second copy
        // Avoid reading unallocated memory
        if(chr>26)
        {
          printf("Warning: Does not know how to color label %u\n", L[kk]);
          chr = 26;
        }
        r = (double) cmap[3*chr]/255.0;
        g = (double) cmap[3*chr+1]/255.0;
        b = (double) cmap[3*chr+2]/255.0;
    }

    sprintf(line, "<marker id=\"%zu\" x=\"%.3f\" y=\"%.3f\" z=\"%.3f\" r=\"%f\" g=\"%f\" b=\"%f\" radius=\"%f\" />\n", 
        kk,
        D[kk*3], D[kk*3+1], D[kk*3+2],
        r, g, b,
        radius);
    gzwrite(zf, line, strlen(line));
  }

  for(size_t kk = 0; kk<NP; kk++)
  {

    size_t A = P[kk*2];
    size_t B = P[kk*2+1];

    double r = .5;
    double g = .5;
    double b = .5;

    if(A==B+1 || A+1==B)
    { // adjacent
      r = .8;
      g = .8;
      b = 0;
    }

    sprintf(line, "<link id1=\"%u\" id2=\"%u\" r=\"%f\" g=\"%f\" b=\"%f\" radius=\"%f\"/>\n",
        P[2*kk], P[2*kk+1], 
        r, g, b,
        radius/3);    
    gzwrite(zf, line, strlen(line));
  }


  sprintf(line, "</marker_set>\n");
  gzwrite(zf, line, strlen(line));

  gzclose(zf);

  free(line);
  return(0);
}




int cmmwrite(char * fname, double * D, size_t nD, double radius, uint32_t * P, size_t NP, uint8_t * L)
{

  FILE * f = fopen(fname, "w");


  if(f == NULL)
  {
    fprintf(stderr, "Unable to open %s\n", fname);
    return -1;
  }

  fprintf(f, "<marker_set name=\"dump\">\n");


  for(size_t kk = 0; kk<nD; kk++)
  {

    double r = 1;
    double g = 0;
    double b = 0;

    if(L != NULL)
    {
        uint8_t chr = L[kk];
        chr = chr % 32; // use 33, 34, ... for the second copy
        // Avoid reading unallocated memory
        if(chr>26)
        {
          printf("Warning: Does not know how to color label %u\n", L[kk]);
          chr = 26;
        }
        r = (double) cmap[3*chr]/255.0;
        g = (double) cmap[3*chr+1]/255.0;
        b = (double) cmap[3*chr+2]/255.0;
    }


    fprintf(f, "<marker id=\"%zu\" x=\"%.3f\" y=\"%.3f\" z=\"%.3f\" r=\"%f\" g=\"%f\" b=\"%f\" radius=\"%f\" />\n", 
        kk,
        D[kk*3], D[kk*3+1], D[kk*3+2],
        r, g, b,
        radius);
  }

  for(size_t kk = 0; kk<NP; kk++)
  {

    size_t A = P[kk*2];
    size_t B = P[kk*2+1];

    double r = .5;
    double g = .5;
    double b = .5;

    if(A==B+1 || A+1==B)
    { // adjacent
      r = .8;
      g = .8;
      b = 0;
    }


    fprintf(f, "<link id1=\"%u\" id2=\"%u\" r=\"%f\" g=\"%f\" b=\"%f\" radius=\"%f\"/>\n",
        P[2*kk], P[2*kk+1], 
        r, g, b,
        radius/3);    
  }

  fprintf(f, "</marker_set>\n");

  fclose(f);
  return(0);
}
