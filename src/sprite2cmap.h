#ifndef _sprite2cmap_h_
#define _sprite2cmap_h_


/* Read .cluster files associated with SPRITE
 *
 * TODO:
 *  - Write more to the log file, i.e., version, date, etc
 *  - Allow creation of heat maps for sub-regions, i.e.
 *    chr1:100,000,000-chr1:120,000,000 vs chr7:10,000,000-chr1:20,000,000
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <libgen.h>
#include <setjmp.h>


/* Command line interface */
int sprite2cmap(int argc, char ** argv);

#endif
