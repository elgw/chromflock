#ifndef _cmmwrite_h_
#define _cmmwrite_h_

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

int cmmwrite(char * fname,
    double * D, // 3xnD list with dot coordinates
    size_t nD, // number of beads/dots
    double radius, // bead radius
    uint32_t * P, // 2xnP list of pairs
    size_t nP, // Number of pairs
    uint8_t * L); // nD long list of labels

// Same as above, but writes a compressed file using zlib
int cmmwritez(char * fname, double * D, size_t nD, double radius, uint32_t * P, size_t nP, uint8_t * L);

#endif
