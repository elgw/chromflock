#pragma once

/**
 * @file mflock.h
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"


#ifdef SDL
#include <pthread.h>
#include "liveview.h"
#endif

#include "cf_version.h"
#include "cf_util.h"
#include "cmmwrite.h"
#include "ellipsoid.h"
#include "functional.h"
#include "wio.h"
#include "contact_pairs_io.h"

/** @breif the command line interface to mflock */
int mflock(int argc, char ** argv);
