#ifndef chromflock_init_h
#define chromflock_init_h

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "oscp.h"

int chromflock_init(void);

#endif
