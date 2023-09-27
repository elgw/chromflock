##
# Purpose:
#    For building chromflock.
#
# NOTES:
# - Please use the makedeb-ubuntu_2204.sh or similar scripts to
#   create a deb file for installation.
#
# WARNING:
#   Do not use --finite-math-only (or --ffast-math or -Ofast) since
#   both aflock and mflock uses isfinite().

CC = gcc -std=gnu99
CFLAGS=-Wall -Wextra
LDFLAGS=

DEBUG?=0

ifeq ($(DEBUG),1)
CFLAGS +=-g3 \
-DNOMATLAB
else
CFLAGS += -O3 \
-DNDEBUG \
-fno-signed-zeros \
-fno-trapping-math \
-fno-math-errno
LDFLAGS += -flto
endif

#
# Inject some information in the binaries
#
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	MANPATH=/usr/local/share/man/man1
endif

CC_VERSION = "$(shell cc --version | head -n 1)"
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"

CFLAGS += -DCC_VERSION=\"$(CC_VERSION)\"
CFLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"

#
# Standard libraries
#

LDFLAGS += -lm  -lpthread  -ldl

# Library: Z
CFLAGS+=`pkg-config zlib --cflags`
LDFLAGS+=`pkg-config zlib --libs`

# Library: LUA
CFLAGS+=-Isrc/lua-5.3.5/src
LDFLAGS+=-Lsrc/lua-5.3.5/src -llua

# Library: Cairo
CFLAGS+=`pkg-config cairo --cflags`
LDFLAGS+=`pkg-config cairo --libs`

# Library: SDL2
SDL?=1
ifeq ($(SDL),1)

CFLAGS += `pkg-config sdl2 --cflags`
CFLAGS += -DSDL
LDFLAGS += `pkg-config sdl2 --libs`

endif

all: bin/chromflock \
     bin/aflock \
     bin/mflock \

chromflock_files=src/chromflock.c \
src/cc2cpm.c \
src/string2any.c \
src/any2string.c \
src/oscp.c \
src/sprite2cmap.c \
obj/chromflock_init.o \
obj/balance.o


## Targets

bin/chromflock: $(chromflock_files)
	$(CC) $(CFLAGS) $(chromflock_files) $(LDFLAGS) -o bin/chromflock

bin/cmmfilter:
	$(CC) $(CFLAGS)  `xml2-config --cflags` src/cmmfilter.c  `xml2-config --libs` $(LDFLAGS) -o bin/cmmfilter

mflock_files = src/mflock.c \
src/functional.c \
src/cmmwrite.c \
src/wio.c \
src/hsvrgb.c \
obj/ellipsoid.o \


bin/mflock: $(mflock_files) makefile
	$(CC) $(CFLAGS) $(mflock_files) -o bin/mflock $(LDFLAGS)

aflock_files = src/aflock.c src/wio.c src/oscp.c obj/ellipsoid.o

bin/aflock: $(aflock_files) makefile
	$(CC) $(CFLAGS) $(aflock_files) -o bin/aflock $(LDFLAGS)

obj/chromflock_init.o: src/chromflock_init.c
	$(CC) -c $(CFLAGS) src/chromflock_init.c -o obj/chromflock_init.o

obj/ellipsoid.o: src/ellipsoid.c
	$(CC) -c $(CFLAGS) src/ellipsoid.c -o obj/ellipsoid.o

obj/balance.o: src/balance.c
	$(CC) -c $(CFLAGS) src/balance.c -o obj/balance.o
