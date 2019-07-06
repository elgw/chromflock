
# Up use a more recent version of gcc than the one provided by Ubuntu, see
# http://tuxamito.com/wiki/index.php/Installing_newer_GCC_versions_in_Ubuntu
#cc = gcc
cc = gcc
#cc = clang

CC_VERSION = "$(shell cc --version | head -n 1)"
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"

f1 = -DCC_VERSION=\"$(CC_VERSION)\"
f2 = $(f1) -DGIT_VERSION=\"$(GIT_VERSION)\"

cflags_dbg = $(f2) -DNOMATLAB -g -Wall

#cflags = -O3 -flto -march=native -Wall -DNOMATLAB -DGSL_RANGE_CHECK_OFF -DNDEBUG -std=gnu99
cflags = $(f2) -O3 -flto -march=native -Wall -DNDEBUG -std=gnu99 -fno-signed-zeros -fno-trapping-math -fno-math-errno -Wextra
# WARNING: Do not use --finite-math-only (or --ffast-math or -Ofast) since
# both aflock and mflock uses isfinite().

ldflags = -lm `pkg-config zlib --libs` -Wall -lpthread

cflagssdl = $(cflags) `pkg-config sdl2 --cflags` 
ldflagssdl = $(ldflags) `pkg-config sdl2 --libs` -lpthread


all: mflock aflock gencm string2any

mflock: 
	$(cc) $(cflags) src/mflock.c src/functional.c src/cmmwrite.c src/wio.c -o bin/mflock $(ldflags) 

mflocksdl: 
	$(cc) $(cflagssdl) src/mflock.c src/functional.c src/cmmwrite.c src/hsvrgb.c src/wio.c -DSDL $(ldflagssdl) -o bin/mflock 

mflock_dbg: 
	$(cc) $(cflags_dbg) src/mflock.c src/functional.c src/cmmwrite.c src/wio.c -o bin/mflock $(ldflags) 

gencm:
	$(cc) $(cflags) src/gen_cm.c -o bin/gen_cm $(ldflags)

aflock:
	$(cc) $(cflags) src/aflock.c src/wio.c -o bin/aflock $(ldflags)

aflock_dbg:
	$(cc) $(cflags_dbg) src/aflock.c src/wio.c -o bin/aflock $(ldflags)

string2any:
	$(cc) $(cflags) src/string2any.c -o bin/string2any $(ldflags)

 
