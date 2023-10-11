#ifndef cf_util_h_
#define cf_util_h_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/** @brief returns the time at the moment
 *
 * @return A string with the current time in the format
 * "YYYY-mm-dd HH:mm:ss" which has to be freed by the caller
 * Example: 2023-10-11 10:42:49
 */
char * cf_timestr();

#endif
