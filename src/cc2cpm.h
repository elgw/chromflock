#ifndef __cc2cpm_h__
#define __cc2cpm_h__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <pthread.h>

#include "balance.h"
#include "cf_version.h"

#define MODE_DEFAULT 0
#define MODE_EQ 1

/* Command line interface */
int cc2cpm(int argc, char ** argv);

#endif
