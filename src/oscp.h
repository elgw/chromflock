#ifndef __oscp_h__
#define __oscp_h__

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <copyfile.h>
#else
#include <sys/sendfile.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

/** @brief Copy a file using operating system calls
 *
 * Unfortunately it is not part of the c
 * standard library and hence we have to make a difference between
 * linux and unix-based systems.
 *
 * TODO Only return fail or not. Also check that all bytes were written.
 * TODO: move to util.c/h
*/

int oscp(const char * source, const char * destination);

#endif
