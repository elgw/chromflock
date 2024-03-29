#pragma once

/**
 * @file cf_util.h
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __linux__
#include <dlfcn.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#endif

/** @brief returns the time at the moment
 *
 * @return A string with the current time in the format
 * "YYYY-mm-dd HH:mm:ss" which has to be freed by the caller
 * Example: 2023-10-11 10:42:49
 */
char * cf_timestr();

/** @brief Limit the memory available to chromflock
 *
 * So that malloc, calloc, etc actually returns NULL at some point
 *
 * @param max_bytes
 * When argument == 0:
 *    Limit the memory to the amount of free memory
 *    right now
 * When argument > 0
 *     Limit memory to the number of bytes specified
 * Only implemented for linux
 * @return EXIT_SUCCESS if the memory limit was set, else EXIT_FAILURE
 */
int limit_mem(size_t max_bytes);


/** @brief delta t
 *
 * Returns the difference in time (s)
 */
double clockdiff(struct timespec* start,
                 struct timespec * finish);

/** @brief get file size in bytes
 * @return -1 on failure;
 */
int64_t cf_file_size(const char * filename);
