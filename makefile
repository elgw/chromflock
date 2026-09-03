##
# Purpose:
#    For building chromflock.
#
# It is recommended that cmake is used (CMakeLists.txt) for most
# users.

CC = gcc -std=gnu11
CFLAGS=-Wall -Wextra -D_FILE_OFFSET_BITS=64
LDFLAGS=

DEBUG?=0
SAN?=0
ANA?=0

ifeq ($(DEBUG),1)
CFLAGS += -g3 \
-DNOMATLAB \
-pedantic
ifeq ($(ANA), 1)
CFLAGS+=-fanalyzer
endif
ifeq ($(SAN), 1)
CFLAGS+=-fsanitize=address
endif
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
CFLAGS += `pkg-config zlib --cflags`
LDFLAGS += `pkg-config zlib --libs`

# Library: LUA
CFLAGS += -Isrc/lua-5.3.5/src
LDFLAGS += -Lsrc/lua-5.3.5/src -llua

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
obj/balance.o \
obj/cf_util.o \
obj/ellipsoid.o \
obj/npio.o \
obj/contact_pairs_io.o \
obj/gzl.o \
obj/con2mflock.o \
obj/collide3_f32.o

## Targets

headers=src/*.h

bin/chromflock: $(chromflock_files) $(headers)
	$(CC) $(CFLAGS) $(chromflock_files) $(LDFLAGS) -o bin/chromflock

bin/cmmfilter:
	$(CC) $(CFLAGS)  `xml2-config --cflags` src/cmmfilter.c  `xml2-config --libs` $(LDFLAGS) -o bin/cmmfilter

mflock_files = src/mflock_cli.c \
obj/ellipsoid.o \
src/mflock.o \
src/functional.o \
src/cmmwrite.o \
src/wio.o \
src/hsvrgb.o \
src/liveview.o \
obj/contact_pairs_io.o \
obj/cf_util.o \
obj/npio.o \
obj/ddict.o \
obj/collide3_f64.o \


bin/mflock: $(mflock_files) makefile
	$(CC) $(CFLAGS) $(mflock_files) -o bin/mflock $(LDFLAGS)

aflock_files = src/aflock.c \
src/wio.c \
src/oscp.c \
obj/ellipsoid.o \
obj/cf_util.o \
obj/contact_pairs_io.o \
obj/npio.o \
obj/collide3_f32.o

bin/aflock: $(aflock_files) makefile
	$(CC) $(CFLAGS) $(aflock_files) -o bin/aflock $(LDFLAGS)

SRCDIR = src
TXTDIR = src/txt
OBJDIR = obj

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(OBJDIR)/%.o : $(SRCDIR)/npio/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<


#TXTHEADERS: FORCE
#	find src/txt -name "*.txt" -execdir xxd -i {} {}.h \;

FORCE:
