#pragma once

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/** @brief Create a chimera CMM file
 *
 * The file can be opened with https://www.cgl.ucsf.edu/chimera/
 *
 * The color map has 27 quite distinct colors and label l will have
 * color_id = l % 32.  if color_id > 26, it will be set to
 * 26. I.e. for diploid structures it is suggested that you label the
 * first copies from 1.... and the second copies from 33 ...
 *
 * Here is an example from
 * https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/volumepathtracer/volumepathtracer.html#markerfiles
 * showing how a chimera file could look like:
 * <marker_set name="marker set 1">
 * <marker id="1" x="-6.1267" y="17.44" z="-3.1338"  radius="0.35217"/>
 * <marker id="2" x="1.5395" y="16.277" z="-3.0339" r="0" g="1" b="1"
 * radius="0.5" note="An example note"/>
 * <link id1="2" id2="1" r="1" g="1" b="0" radius="0.17609"/>
 * </marker_set>
 *
 * @param D A 3xN list with dot coordinates
 * @param N Number of dots
 * @param radius bead radius
 * @param P 2xnP list of connected beads
 * @param nP number of pairs
 * @param L bead labels which will be used for coloring using the built-in colormap
 *
 * All input parameters are required.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 */
int cmmwrite(const char * fname,
             const double * D,
             size_t N,
             double radius,
             const uint32_t * P,
             size_t nP,
             const uint8_t * L);

/* @brief as cmmwrite but will write to a cmm.gz file using libz
 */
int cmmwritez(const char * fname,
              const double * D,
              size_t nD,
              double radius,
              const uint32_t * P,
              size_t nP,
              const uint8_t * L);
