#ifndef any2str_h_
#define any2str_h_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* TODO:
 * If no data is given, i.e., if argc==3,
 * read from stdin
 */

/* Command line interface */
int any2string(int argc, char ** argv);

#endif
