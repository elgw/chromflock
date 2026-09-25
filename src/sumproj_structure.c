#include "sumproj_structure.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cf_util.h"

#include "txt/sumproj_help.txt.h"

static void
blit_gaussian(float * X,
              int imw,
              float fx, float fy,
              float sigma)
{

    int w = round(2*sigma + 1);
    int xmin = fmax(0, round(fx) - w);
    int xmax = fmin(round(fx) + w, imw-1);
    int ymin = fmax(0, round(fy) - w);
    int ymax = fmin(round(fy) + w, imw -1);
    for(int yy = ymin; yy <= ymax; yy++){
        float ey0 = erf( ((float) yy - fy - 0.5)/sqrt(2)/sigma );
        float ey1 = erf( ((float) yy - fy + 0.5)/sqrt(2)/sigma );
        for(int xx = xmin; xx <= xmax; xx++){
            // integrate over the pixel
            float ex0 = erf( ((float) xx - fx - 0.5)/sqrt(2)/sigma );
            float ex1 = erf( ((float) xx - fx + 0.5)/sqrt(2)/sigma );
            X[xx + yy*imw] += ex0*ey0 + ex1*ey1 - ex0*ey1 - ex1*ey0;
            // The lazy version would be
            //float d = sqrt(   pow((float) fx - xx, 2)
            //+ pow((float) fy - yy, 2));
            //X[xx + yy*imw] += expf( - d*d ); // vs central value

        }
    }
}

// TODO would be nice to expose via Python interface
float *
sumproj_coords(float * X, // table of point coordinates
               size_t N, // number of rows
               size_t row_stride, // elements per row
               size_t imw, // side length of output image
               double sigma, // size of spots
               int verbose)
{


    // Figure out bounding box
    float min[2];
    min[0] = X[0]; min[1] = X[1];
    float max[2];
    max[0] = X[0]; max[1] = X[1];

    for(size_t kk = 0; kk < N; kk++)
    {
        float * P = X + kk*row_stride;
        min[0] = fmin(min[0], P[0]);
        min[1] = fmin(min[1], P[1]);
        max[0] = fmax(max[0], P[0]);
        max[1] = fmax(max[1], P[1]);
    }
    if(verbose > 2){
        printf("BBX: [%f, %f] x [%f, %f]\n", min[0], max[0], min[1], max[1]);
    }
    float center[2];
    center[0] = (min[0] + max[0])/2.0;
    center[1] = (min[1] + max[1])/2.0;

    // Figure out bounding box

    float * R = calloc(imw*imw, sizeof(float));
    if(R == NULL){ return NULL; }

    for(size_t kk = 0; kk < N; kk++)
    {
        // transformation to image coordinates
        float * P = X + kk*row_stride;
        // (fx, fy) sub pixel location of the spot center
        float fx = (P[0] - min[0] - center[0]) / (max[0]-min[0]) * ((double) imw - 1.0);
        float fy = (P[1] - min[1] - center[1]) / (max[1]-min[1]) * ((double) imw - 1.0);

        blit_gaussian(R, imw, fx, fy, sigma);
    }


    return R;
    // done :)
}

static int
load_structure_npy(const char * infile, float ** X, size_t * n, size_t * stride, int verbose)
{
    if(verbose > 2){
        printf("Trying to load data from %s\n", infile);
    }
    npio_t *  coords = npio_load(infile);
    if(coords == NULL){
        fprintf(stderr, "Could not load %s as an npy file\n", infile);
    }

    if(coords->ndim != 2){
        fprintf(stderr, "The list of points must be a 2D array\n");
        goto fail1;
    }

    if(coords->dtype != NPIO_F32){
        fprintf(stderr, "The list of points be 32-bit floats\n");
        goto fail1;
    }

    X[0] = (float*) coords->data;
    coords->data = NULL;
    *n = coords->shape[0];
    *stride = coords->shape[1];
    npio_free(coords);
    return 0;

 fail1:
    npio_free(coords);
    return -1;
}

static int
load_structure_3dg(const char * infile, float ** _X, size_t * n, size_t * stride, int verbose)
{
    // Example line:
    // chr10_0 40000   -99.5743        21.6467 -10.5219
    size_t n_lines = 0;
    FILE * fid = fopen(infile, "r");
    if(fid == NULL){
        fprintf(stderr, "Unable to open %s\n", infile);
        return -1;
    }
    char * linep = NULL;
    size_t linesize = 0;
    while(getline(&linep, &linesize, fid) >= 0){
        n_lines++;
    }
    rewind(fid);
    if(verbose > 1){
        printf("found %zu lines\n", n_lines);
    }

    if(n_lines < 1){
        free(linep);
        return -1;
    }

    float * X = calloc(n_lines*3, sizeof(float));
    if(X == NULL){
        free(linep);
        return -1;
    }
    _X[0] = X;
    size_t npoint = 0;
    while(getline(&linep, &linesize, fid) >= 0){
        char *token = strtok(linep, "\t ");
        size_t nt = 0;
        assert(npoint < n_lines);
        while( token ){
            nt++;
            if((nt >= 3) && (nt <= 5)){
                X[3*npoint + nt - 3] = atof(token);
            }
            // No proper error checking here, but if we found
            // the expected number of white space separated columns
            // it is a go.
            if(nt == 5){
                npoint++;
            }
            token = strtok(NULL, "\t ");
        }
    }

    *n = npoint;
    *stride = 3;

    free(linep);
    fclose(fid);
    return 0;
}

static int
load_structure(const char * infile, float ** X, size_t * n, size_t * stride, int verbose)
{
    if(npy_extension(infile)){
        return load_structure_npy(infile, X, n, stride, verbose);
    } else {
        return load_structure_3dg(infile, X, n, stride, verbose);
    }
}

typedef struct {
    u32 imsize;
    float sigma;
    char **files;
    int nfiles;
    int verbose;
} config;

config * config_new(int argc, char ** argv){
    config * conf = calloc(1, sizeof(config));
    conf->sigma = 1.5;
    conf->imsize = 400;
    conf->verbose = 1;

    struct option options[] = {
        {"sigma", required_argument,  NULL, 's'},
        {"size",        required_argument,  NULL, 'w'},
        {"help", no_argument, NULL, 'h'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};

    int ch;
    while((ch = getopt_long(argc, argv, "s:w:hv", options, NULL)) != -1)
    {
        switch(ch) {
        case 's':
            conf->sigma = atof(optarg);
            break;
        case 'w':
            conf->imsize = atol(optarg);
            break;
        case 'h':
            printf("%s", sumproj_help_txt);
            free(conf);
            return NULL;
        case 'v':
            conf->verbose = atoi(optarg);
        }
    }
    conf->files = argv + optind;
    conf->nfiles = argc - optind;
    return conf;
}

int sumproj_structure(int argc, char ** argv)
{
    config * conf = config_new(argc, argv);
    if(conf == NULL){
        return -1;
    }

    if(conf->nfiles == 0){
        fprintf(stderr, "sumproj_structure: no input file given\n");
        goto fail1;
    }

    for(int kk = 0; kk < conf->nfiles; kk++){
        const char * infile = conf->files[kk];

        float * X = NULL;
        size_t n = 0;
        size_t stride = 0;

        if(load_structure(infile, &X, &n, &stride, conf->verbose)){
            printf("Unable to load %s\n", infile);
            goto fail1;
        }

        char * outfile = calloc(32 + strlen(infile), 1);
        if(outfile == NULL){
            fprintf(stderr, "ERROR: sumproj_structure, L%d\n", __LINE__);
            free(X);
            goto fail1;
        }
        sprintf(outfile, "%s_sumz.npy", infile);
        if(conf->verbose > 1){
            printf("%s -> %s\n", infile, outfile);
        }
        float * raster = sumproj_coords(X, n, stride,
                                        conf->imsize, conf->sigma,
                                        conf->verbose);
        free(X);
        if(raster == NULL){
            free(outfile);
            goto fail1;
        }

        int out_shape[2] = {conf->imsize, conf->imsize};
        npio_write(outfile, 2, out_shape, (void*) raster, NPIO_F32, NPIO_F32);
        free(raster);
        free(outfile);
    }

    free(conf);
    return 0;

 fail1:
    free(conf);
    return -1;
}
