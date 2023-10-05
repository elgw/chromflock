CC=gcc -std=gnu99 -Wall -Wextra

CFLAGS = -DNOMATLAB -DMAIN

ifeq ($(DEBUG),1)
CFLAGS += -g3 \
-DNOMATLAB \
-fanalyzer
else
CFLAGS += -O3 -flto
endif

CFLAGS+=-DCC_VERSION='"..."'
CFLAGS+=-DGIT_VERSION='"..."'

CFLAGSdbg = -g3 -Wall -DNOMATLAB -DMAIN


CFLAGSlv = -std=gnu99 -march=native -Wall -g `pkg-config sdl2 --cflags`
LDFLAGSlv = -lm -lpthread `pkg-config sdl2 --libs` -Wall -lpthread

CFLAGScmm = $(CFLAGS) `pkg-config cairo --cflags`
LDFLAGScmm += -lz `pkg-config --libs cairo`

LDFLAGS=-lm

all: functional_ut cmmwrite liveview hsvrgb ellipsoid functional_optimal_nDiv wio_ut

wio_ut: wio.c wio.h makefile
	$(CC) $(CFLAGS) wio_ut.c wio.c -o wio_ut -lz

ellipsoid:
	$(CC) $(CFLAGS) ellipsoid.c -lm -o ellipsoid_ut -Dstandalone

functional_ut:
	$(CC) $(CFLAGS) ellipsoid.c functional_ut.c functional.c -lm -o functional_ut

functional_optimal_nDiv:
	$(CC) $(CFLAGS) functional_optimal_nDiv.c functional.c ellipsoid.c  -lm -o functional_optimal_nDiv

cmmwrite:
	$(CC) cmmwrite_ut.c  $(CFLAGScmm) $(LDFLAGScmm) -o cmmwrite
	# ./cmmwrite
	# gunzip testz.cmm.gz
	# diff test.cmm testz.cmm

FILESlv = liveview.c liveview_ut.c ellipsoid.c hsvrgb.c
liveview:
	$(CC) $(CFLAGSlv) $(FILESlv) $(LDFLAGSlv) -o liveview

hsvrgb:
	$(CC) $(CFLAGS) hsvrgb_ut.c -lm -o hsvrgb


FILES=cc2cpm.c cc2cpm_cli.c balance.c
cc2cpm_cli: $(FILES)
	$(CC) $(CFLAGS) $(FILES) $(LDFLAGS) -o cc2cpm_cli
