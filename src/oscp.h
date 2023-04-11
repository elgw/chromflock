#ifndef __oscp_h__
#define __oscp_h__

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

/* Copy file using os calls Unfortunately it is not part of the c
 * standard library and hence we have to make a difference between
 * linux and unix-based systems.
 * Return value:
 * - On linux: see sendfile
 * - On MacOS: see fcopyfile
 * TODO Only return fail or not. Also check that all bytes were written.
 * TODO: move to util.c/h
*/

int oscp(const char * , const char * );

#endif
