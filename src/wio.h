#ifndef wio_h
#define wio_h

/* Read/Write contact indication matrices stored as uint8 */

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include <assert.h>

/* Read an array stored as uint8_t. Use lz if the file ends with '.gz'
 * returns NULL on failure
 * WARNING: uses the last 4 bytes of 'gz' files to determine file size
 * that fail if the file is  2^32 bytes or larger uncompressed.
 */
void * wio_read(char * fileName, size_t * nel);

/* Write an array stored as uint8_t. Use lz if the file ends with '.gz'
 * returns 0 on sucess
 */
int wio_write(char * fileName, size_t nel, void * data);

/* Unit tests called by wio_ut.c
 */
int wio_ut(int argc, char ** argv);

#endif
