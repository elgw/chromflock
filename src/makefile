CC=gcc

CFLAGS = -std=gnu99 -flto -O3 -march=native -mtune=native -Wall -DNOMATLAB -DMAIN
CFLAGS+=-DCC_VERSION='"..."'
CFLAGS+=-DGIT_VERSION='"..."'

CFLAGSdbg = -g3 -Wall -DNOMATLAB -DMAIN
CFLAGScmm = $(CFLAGS) `pkg-config cairo --CFLAGS --libs`

CFLAGS_lv = -std=gnu99 -march=native -Wall -g `pkg-config sdl2 --CFLAGS`
ldflags_lv = -lm -lpthread `pkg-config sdl2 --libs` -Wall -lpthread

LDFLAGS=-lm

all: functional_ut cmmwrite liveview hsvrgb ellipsoid hic2a functional_optimal_nDiv wio_ut

wio_ut: wio.c wio.h makefile
	$(CC) $(CFLAGS) wio_ut.c wio.c -o wio_ut -lz

ellipsoid:
	$(CC) $(CFLAGS) ellipsoid.c -lm -o ellipsoid_ut -Dstandalone

functional_ut:
	$(CC) $(CFLAGS) ellipsoid.c functional_ut.c functional.c -lm -o functional_ut

functional_optimal_nDiv:
	$(CC) $(CFLAGS) functional_optimal_nDiv.c functional.c ellipsoid.c  -lm -o functional_optimal_nDiv

cmmwrite:
	$(CC) cmmwrite_ut.c  $(CFLAGScmm) -lz -o cmmwrite
	# ./cmmwrite
	# gunzip testz.cmm.gz
	# diff test.cmm testz.cmm

liveview:
	$(CC) $(CFLAGS_lv) liveview_ut.c liveview.c $(ldflags_lv) -o liveview

hsvrgb:
	$(CC) $(CFLAGS) hsvrgb_ut.c -lm -o hsvrgb

mtwist:
	$(CC) $(CFLAGS) -c mtwist.c
	$(CC) $(CFLAGS) mtwist_ut.c mtwist.o -o mtwist

FILES=cc2cpm.c cc2cpm_cli.c balance.c
cc2cpm_cli: $(FILES)
	$(CC) $(CFLAGS) $(FILES) $(LDFLAGS) -o cc2cpm_cli
