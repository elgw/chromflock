#include "cmmwrite.h"

typedef uint8_t u8;

static const u8 cmap0[] =
{ 255,255,255,
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

// Default colormap. 1 ... 24 colored as well as 33 .. 45, everything
// else is white

static u8 * default_cmap(void)
{
    u8 * cmap = calloc(256*3, sizeof(u8));
    for(int kk = 0; kk < 256*3; kk++)
    {
        cmap[kk] = 255;
    }
    for(int kk = 0; kk < 27; kk++)
    {
        memcpy(cmap + 3*kk, cmap0 + 3*kk, 3);
        memcpy(cmap + 3*(kk+32), cmap0 + 3*kk, 3);
    }
    return cmap;
}

// if a colormap is provided,
// it should have 256 colors (set whatever isn't used to whatever you want)
static int
cmmwrite_general(const char * fname,
                 const double * D, size_t nD, double radius,
                 const uint32_t * P, size_t NP,
                 const uint8_t * L,
                 bool useGZ,
                 const u8 * cmap)
{
    u8 * cmap_default = NULL;
    if(cmap == NULL)
    {
        cmap_default = default_cmap();
        cmap = cmap_default;
    }
    if(L == NULL)
    {
        fprintf(stderr, "ERROR: cmmwrite needs a non-null list of labels\n");
        exit(EXIT_FAILURE);
    }

    gzFile zf = NULL;
    FILE * f = NULL;;

    if(useGZ)
    {
        zf = gzopen(fname, "wb");
        if(zf == Z_NULL)
        {
            fprintf(stderr, "Unable to open %s\n", fname);
            return EXIT_FAILURE;
        }
    } else {
        f = fopen(fname, "w");
        if(f == NULL)
        {
            fprintf(stderr, "cmmwrite: Unable to open %s for writing\n", fname);
            return EXIT_FAILURE;
        }
    }

    char * line = malloc(1024*sizeof(char));
    if(line == NULL)
    {
        if(useGZ)
        {
            gzclose(zf);
        } else
        {
            fclose(f);
        }
        return EXIT_FAILURE;
    }

    sprintf(line, "<marker_set name=\"dump\">\n");
    if(useGZ)
    {
        gzwrite(zf, line, strlen(line));
    } else {
        fprintf(f, "%s", line);
    }


    for(size_t kk = 0; kk<nD; kk++)
    {
        double r = 1;
        double g = 0;
        double b = 0;

        uint8_t chr = L[kk];

        r = (double) cmap[3*chr]/255.0;
        g = (double) cmap[3*chr+1]/255.0;
        b = (double) cmap[3*chr+2]/255.0;

        //printf("%u(%d) -> %f, %f, %f\n", chr, chr % 32, r, g, b);

        sprintf(line, "<marker id=\"%zu\" x=\"%.3f\" y=\"%.3f\" z=\"%.3f\" r=\"%f\" g=\"%f\" b=\"%f\" radius=\"%f\" />\n",
                kk,
                D[kk*3], D[kk*3+1], D[kk*3+2],
                r, g, b,
                radius);
        if(useGZ)
        {
            gzwrite(zf, line, strlen(line));
        } else {
            fprintf(f, "%s", line);
        }
    }

    /* Write links between beads */
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
        if(useGZ)
        {
            gzwrite(zf, line, strlen(line));
        } else {
            fprintf(f, "%s", line);
        }
    }


    sprintf(line, "</marker_set>\n");
    if(useGZ)
    {
        gzwrite(zf, line, strlen(line));
    } else {
        fprintf(f, "%s", line);
    }

    if(useGZ)
    {
        gzclose(zf);
    } else {
        fclose(f);
    }

    free(line);


    free(cmap_default);

    return EXIT_SUCCESS;
}


int cmmwritez(const char * fname,
              const double * D, size_t nD, double radius,
              const uint32_t * P, size_t NP,
              const uint8_t * L,
              const uint8_t * cmap)
{
    return cmmwrite_general(fname, D, nD, radius, P, NP, L, 1, cmap);
}


int cmmwrite(const char * fname,
             const double * D, size_t nD, double radius,
             const uint32_t * P, size_t NP,
             const uint8_t * L,
             const uint8_t * cmap)
{
    return cmmwrite_general(fname, D, nD, radius, P, NP, L, 0, cmap);
}
