#pragma once

/* Simple line-by line reader for zlib compressed files
 * (.gz). Developed against zlib.h version 1.3.
 *
 * Main Caveat/Feature:
 *
 * - Works with non-compressed files as well.
 *
 * - It will fail (gracefully) if it find a line larger than the
 *   max_line_size.
 *
 * Version 1.0.0
 * Erik Wernersson 2026-01-13
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <zlib.h>

typedef struct gzl__state gzl_state;

/* Open a file for reading lines at most max_line_size long
 *
 * On error: returns NULL
 */

gzl_state *
gzl_open(const char * filename, int64_t max_line_size);

/*
 *
 */

void
gzl_destroy(gzl_state * gzl);


/*
 * Returns the next line in the file or NULL. A returned string is always
 * terminated by '\0' (and will not contain any newline character).
 *
 * The returned string is owned by the gzl object and should not be
 * freed. It is safe to modify the strlen bytes of the returned
 * string. The returned string can only be used until the next call
 * to gzl_get_line.
 *
 * If the return value is NULL, please check the error value:
 *
 *   0 Normal exit, i.e. end of file was reached
 *   1 Found a line longer than the max_line_size
 *   2 internal library error, something unexpected happened, please
 *     file a bug report.
 * < 0 indicates a zlib error code.
 *
 * Corresponding, human readable error message can be retrived by
 *     gzl_get_error_string.
 *
 * After an error has occured the object can not read anything more
 * from the file and should be destroyed.
 */

const char *
gzl_get_line(gzl_state * gzl, int * gzl_error);

/* Human readable string based on the error state */
const char *
gzl_get_error_string(gzl_state * gzl);

#ifdef __cplusplus
}
#endif
