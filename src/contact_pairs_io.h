#ifndef contact_pairs_h_
#define contact_pairs_h_

/**
 * @file contact_pairs_io.h
 * @brief Read/write contact pairs as raw uint32 files.
 * @author Erik Wernersson
 * @date 2023
 *
 * A contact pair (a,b), where a, and b, are bead numbers, starting
 * from 0, indicates that two beads are in contact. Chromflock
 * requires that a<b for each pair. It further assumes that
 * a_i < a_j, i.e. that the pairs are sorted.
 *
 */

#include <assert.h>
#include <fcntl.h>
#include <math.h>
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

/** @brief Extract contact pairs from a NxN=n_elements matrix
 *         and write them to disk.
 *
 * As a replacement for wio_write
 *
 * @param filename, will be used as in contact_pairs_write.
 * @param n_elements the number of elements in the matrix W
 * @param W a square matrix.
*/
int contact_pairs_write_from_matrix(const char * filename,
                                    uint64_t n_elements ,
                                    const uint8_t * W);

/** @brief Contact indication matrix to list of contacts
 *
 * Note that the diagonal is ignored.
 *
 * @param W a N x N matrix assumed to be symmetric
 * @param N the size of W
 * @param[out] nPairs will be set to the number of pairs found.
 * @return A list of pairs in contact.
 */
uint32_t *
    contact_pairs_from_matrix(const uint8_t * restrict W,
                              const uint64_t N,
                              uint64_t * nPairs);

/** @brief Convert a list of contact pairs to a indication matrix
 *
 * @param CP The contacts
 * @param nPairs is the number of pairs
 * @param N is the number of beads
 * @returns a symmetric NxN matrix where contacts are indicated by a 1 or NULL on failure.
*/
uint8_t * contact_pairs_to_matrix(const uint32_t * CP, uint64_t nPairs, uint64_t N);

/* Unit tests */
int contact_pairs_io_ut(int argc, char ** argv);

#endif
