#ifndef contact_pairs_h_
#define contact_pairs_h_

/**
 * @file contact_pairs_io.h
 * @brief Read/write contact pairs as raw uint32 files.
 * @author Erik Wernersson
 * @date 2023
 */

#include <assert.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <zlib.h>
#include <unistd.h>

/** @brief Read contact pairs from filename
if the file ends with gz, libzg will be used to read it.

@param[in] filename file to read from.
@param[out] nPairs number of contact pairs in the file.
@return NULL on failure, in that case nPairs is undefined.
*/
uint32_t * contact_pairs_read(const char * filename,
                              uint64_t * nPairs);

/** @brief Write contact pairs to filename
@param[in] filename file to write to
@param[in] P contact pairs
@param[in] nPairs number of contact pairs
@return EXIT_FAILURE on failure

if the filename ends with gz, libgz will be used to compress it
*/
int contact_pairs_write(const char * filename,
                        const uint32_t * P,
                        uint64_t nPairs);

/* Unit tests */
int contact_pairs_io_ut(int argc, char ** argv);

#endif
