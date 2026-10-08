#pragma once

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

// Create a chimera CMM file
//
// The file can be opened with https://www.cgl.ucsf.edu/chimera/
//
// The color map has 27 quite distinct colors and label l will have
// color_id = l % 32.  if color_id > 26, it will be set to
// 26. I.e. for diploid structures it is suggested that you label the
// first copies from 1.... and the second copies from 33 ...
//
// Here is an example from
// https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/volumepathtracer/volumepathtracer.html#markerfiles
// showing how a chimera file could look like:
// <marker_set name="marker set 1">
// <marker id="1" x="-6.1267" y="17.44" z="-3.1338"  radius="0.35217"/>
// <marker id="2" x="1.5395" y="16.277" z="-3.0339" r="0" g="1" b="1"
// radius="0.5" note="An example note"/>
// <link id1="2" id2="1" r="1" g="1" b="0" radius="0.17609"/>
// </marker_set>
//
// Parameters:
//  fname output file name. If it ends with .gz the file will be written
//  with libgz
//  D A 3xN list with dot coordinates
//  N Number of dots
//  radius bead radius
//  P 2xnP list of connected beads
//  nP number of pairs
//  L bead labels which will be used for coloring using the built-in colormap
//  cmap A 3 x 256 array containing a RGB colormap. Can be null.
//
// All input parameters are required.
//
// returns EXIT_SUCCESS or EXIT_FAILURE

int cmmwrite(const char * fname,
             const double * D,
             size_t N,
             double radius,
             const uint32_t * P, // backbone contacts
             size_t nP,

             // Contacts
             // If CTI is NULL all contacts are considered enabled
             const uint32_t * CT, // 2 u32 per contact
             const uint32_t n_CT, // number of contacts
             const uint8_t * CTI, // enabled or not indicator

             const uint8_t * L,
             const uint8_t * cmap);
