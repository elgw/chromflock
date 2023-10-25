#pragma once

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

/* Command line interface */
int cc2cpm(int argc, char ** argv);
