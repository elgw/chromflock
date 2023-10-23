#include "cmmwrite.h"

static uint8_t cmap[] = {255,255,255,
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


static int
cmmwrite_general(const char * fname,
                 const double * D, size_t nD, double radius,
                 const uint32_t * P, size_t NP,
                 const uint8_t * L,
                 bool useGZ)
{
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

    int labelWarning = 0;
    for(size_t kk = 0; kk<nD; kk++)
    {
        double r = 1;
        double g = 0;
        double b = 0;

        uint8_t chr = L[kk];
        chr = chr % 32; // use 33, 34, ... for the second copy
        // Avoid reading unallocated memory
        if(chr>26)
        {
            labelWarning = 1;
            chr = 26;
        }
        r = (double) cmap[3*chr]/255.0;
        g = (double) cmap[3*chr+1]/255.0;
        b = (double) cmap[3*chr+2]/255.0;

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

    if(labelWarning)
    {
        fprintf(stderr,
                "Warning: cmmwrite got labels, l, where l %% 32 > 26 "
                "that was not expected and indicates that the function was used "
                "in a non intended way. Please file a bug report!\n");
    }

    return EXIT_SUCCESS;
}


int cmmwritez(const char * fname,
              const double * D, size_t nD, double radius,
              const uint32_t * P, size_t NP,
              const uint8_t * L)
{
    return cmmwrite_general(fname, D, nD, radius, P, NP, L, 1);
}


int cmmwrite(const char * fname,
             const double * D, size_t nD, double radius,
             const uint32_t * P, size_t NP,
             const uint8_t * L)
{
    return cmmwrite_general(fname, D, nD, radius, P, NP, L, 0);
}
