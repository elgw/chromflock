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
#include <stdint.h>

#ifdef __linux__
#include <dlfcn.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#endif

#include "ellipsoid.h"

typedef uint32_t u32;
typedef uint8_t u8;

/* For holding absolute bead position information used by mflock (--absolute) */
typedef struct {
    u32 bead_id;
    float x;
    float y;
    float z;
} bpos;

void bpos_print(FILE * fid, bpos *);

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

/* Returns 1 if the string ends with .npy */
int npy_extension(const char * name);

/* Write an array of bead coordinates to a csv file
 * For each bead, x, y, z and r will be written.
 * The reason for writing the radius is a convenience
 * when the geometry is non-spherical (ellipsoidal)
 * E should be set to NULL when a spherical geometry is used
 *
 * returns 0 on success
 */

int
write_bead_coordinates_to_csv(const char * fname,
                              const double * X,
                              const int64_t nbead,
                              const elli * geometry);

int
write_bead_coordinates_to_npy(const char * fname,
                                  const double * X,
                                  const int64_t nbead,
                                  const elli * geometry);

/* Read nbead rows from a csv
 * Does not expect a header rows
 * Three values are read from row, any extra values are ignored
 * Values are interpreted as x, y, z coordinates of a bead
 *
 * Writes values to either X32 or X64, i.e. exactly one of them should
 * be non-NULL.
 */

int
load_bead_coordinates_from_csv(const char * fname,
                               float * X32, double * X64,
                               const int64_t nbead);

int
load_bead_coordinates_from_npy(const char * fname,
                               float * X32, double * X64,
                               const int64_t nbead);


// Load absolute positions for (certain) beads
//
// On success: sets the number of constraints to nconstraint
//
// On failure: Returns NULL

bpos *
load_bead_apos_from_npy(const char * fname,
                        int * nconstraint);

u8 *
load_bead_labels_from_npy(const char * fname, int * nbead);

u32 *
load_bead_contacts_from_npy(const char * fname, int * ncont);
