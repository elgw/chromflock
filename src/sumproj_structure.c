#include "sumproj_structure.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cf_util.h"
#include "npio.h"

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
               double sigma) // size of spots
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
    printf("BBX: [%f, %f] x [%f, %f]\n", min[0], max[0], min[1], max[1]);
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
load_structure_npy(const char * infile, float ** X, size_t * n, size_t * stride)
{
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
load_structure_3dg(const char * infile, float ** _X, size_t * n, size_t * stride)
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
    printf("%zu lines\n", n_lines);

    float * X = malloc(n_lines*3*sizeof(float));
    if(X == NULL){
        free(linep);
        return -1;
    }
    _X[0] = X;
    size_t npoint = 0;
    while(getline(&linep, &linesize, fid) >= 0){
        char *token = strtok(linep, "\t ");
        size_t nt = 0;
        while( token ){
            nt++;
            if((nt >= 3) && (nt <= 5)){
                X[3*npoint + nt -3] = atof(token);
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
load_structure(const char * infile, float ** X, size_t * n, size_t * stride)
{
    if(npy_extension(infile)){
            return load_structure_npy(infile, X, n, stride);
        } else {
        return load_structure_3dg(infile, X, n, stride);
    }
}

// For one file, either the coordinates from mcflock, i.e.
// coords.npy or from a 3dg file
// Create a 2D image with the sum projection over z
// TODO ? : add custom orientation as input
// TODO ? : select chromosome
int sumproj_structure(int argc, char ** argv)
{
    size_t imw = 400;
    float sigma = 1.0;

    if(argc < 2){
        fprintf(stderr, "sumproj_structure: no input file given\n");
        return 1;
    }

    const char * infile = argv[1];


    if(argc > 2) {
        sigma = atof(argv[2]);
    }
    if(argc > 3){
        imw = atol(argv[3]);
    }

    float * X = NULL;
    size_t n = 0;
    size_t stride = 0;

    if(load_structure(infile, &X, &n, &stride)){
        printf("Unable to load %s\n", infile);
        return -1;
    }

    char * outfile = calloc(32 + strlen(infile), 1);
    if(outfile == NULL){
        fprintf(stderr, "ERROR: sumproj_structure, L%d\n", __LINE__);
        return -1;
    }
    sprintf(outfile, "%s_maxz.npy", infile);

    printf("sigma = %f\n", sigma);
    printf("%s -> %s\n", infile, outfile);


    printf("Projecting %zu points\n", n);


    float * raster = sumproj_coords(X, n, stride,
                                    imw, sigma);
    if(raster == NULL){
        goto fail1;
    }

    int out_shape[2] = {imw, imw};
    npio_write(outfile, 2, out_shape, (void*) raster, NPIO_F32, NPIO_F32);
    free(raster);
    free(outfile);
    return 0;

 fail1:
    free(outfile);
    return -1;

}
