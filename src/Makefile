#
#  File        : Makefile
#                ( Makefile for GNU 'make' utility )
#
#  Description : This Makefile exposes the following entries:
#
#                 . 'all':
#                   Equivalent to 'lib','cli','gimp','gmic_qt','libc','zart' (see below).
#
#                 . 'lib':
#                   C++ API for the G'MIC library. Generate files 'src/libgmic.*'.
#
#                 . 'libc' and 'libc_static':
#                   C API for the G'MIC library. Generate files 'src/libcgmic.*'.
#                   'libc_static' tries to embed most libraries as static.
#
#                 . 'cli', 'cli_shared' and 'cli_static':
#                    G'MIC command line tool. Generate file 'src/gmic[.exe]'
#                    'cli_shared' generates an executable that is dynamically linked with 'libgmic'.
#                    'cli_static' generates an executable where most libraries are statically linked.
#
#                 . 'gimp' and 'gimp_shared':
#                   G'MIC-Qt plug-in for GIMP. Generate file 'gmic-qt/gmic_gimp_qt'.
#                   'gimp_shared' generates a binary that is dynamically linked with 'libgmic'.
#
#                 . 'gmic_qt' and 'gmic_qt_shared':
#                   G'MIC-Qt stand-alone application. Generate file 'gmic-qt/gmic_qt'.
#                   'gimp_shared' generates a binary that is dynamically linked with 'libgmic'.
#
#                 . 'zart':
#                   ZArt interface for real-time processing of videos coming from webcams or files.
#                   Generate file 'zart/zart'.
#
#                ( https://gmic.eu )
#
#  Copyright   : David Tschumperlé
#                ( https://tschumperle.users.greyc.fr/ )
#
#  Licenses    : This file is 'dual-licensed', you have to choose one
#      	          of the two licenses below to apply.
#
#                CeCILL-C
#                The CeCILL-C license is close to the GNU LGPL.
#                ( http://cecill.info/licences/Licence_CeCILL-C_V1-en.html )
#
#            or  CeCILL v2.1
#                The CeCILL license is compatible with the GNU GPL.
#                ( http://cecill.info/licences/Licence_CeCILL_V2.1-en.html )
#
#  This software is governed either by the CeCILL or the CeCILL-C license
#  under French law and abiding by the rules of distribution of free software.
#  You can  use, modify and or redistribute the software under the terms of
#  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA
#  at the following URL: "http://cecill.info".
#
#  As a counterpart to the access to the source code and  rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty  and the software's author,  the holder of the
#  economic rights,  and the successive licensors  have only  limited
#  liability.
#
#  In this respect, the user's attention is drawn to the risks associated
#  with loading,  using,  modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean  that it is complicated to manipulate,  and  that  also
#  therefore means  that it is reserved for developers  and  experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or
#  data to be ensured and,  more generally, to use and operate it in the
#  same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.
#

#---------------------------
# Set OS-specific variables.
#---------------------------

# Try to auto-detect target OS.
OS = $(shell uname)
PREFIX ?= /usr
LIB ?= lib
BIN ?= bin
INCLUDE ?= include
PIC = -fPIC
EXE =
WGET = wget --no-check-certificate --quiet -O
PLUGINDIR = $(shell gimptool-2.0 --gimpplugindir)/plug-ins
BASH_COMPLETIONSDIR = $(shell pkg-config --variable completionsdir bash-completion || echo $(PREFIX)/share/bash-completion/completions)
VERSION = $(shell grep 'gmic_version\ ' gmic.h | tail -c4 | head -c3)
VERSION1 = $(shell grep 'gmic_version\ ' gmic.h | tail -c4 | head -c1)
VERSION2 = $(shell grep 'gmic_version\ ' gmic.h | tail -c3 | head -c1)
VERSION3 = $(shell grep 'gmic_version\ ' gmic.h | tail -c2 | head -c1)
SVERSION=$(VERSION1).$(VERSION2).$(VERSION3)
SRC_PATH=$(shell pwd)

# Check that versions of files 'CImg.h' and 'gmic.h' match.
ifneq (,$(wildcard CImg.h))
CIMG_VERSION = $(shell grep 'cimg_version\ ' CImg.h | tail -c4 | head -c3)
else
CIMG_VERSION = $(VERSION)
endif

ifneq ($(VERSION),$(CIMG_VERSION))
REMOVE_CIMG_H := $(shell rm -f CImg.h)
endif

ifeq ($(OS),Linux)
OS = Unix
endif

ifeq ($(OS),GNU/kFreeBSD)
OS = Unix
endif

ifeq ($(OS),GNU)
OS = Unix
endif

ifeq ($(OS),FreeBSD)
OS = Unix
PREFIX = /usr/local
endif

ifeq ($(OS),DragonFly)
OS = Unix
PREFIX = /usr/local
endif

ifeq ($(OS),NetBSD)
OS = Unix
endif

ifneq (,$(findstring MINGW,$(OS)))
OS = Windows
PREFIX = /mingw64
endif

ifeq ($(OS),Darwin)
ifeq (,$(wildcard /opt/local))
PREFIX = /usr/local
PLUGINDIR = $(HOME)/Library/Application\ Support/Gimp/plug-ins
else
PREFIX = /opt/local
endif
WGET = curl -k -L --silent -o
endif

ifeq ($(OS),Windows)
EXE = .exe
PIC =
endif

ifneq (,$(findstring g++,$(CXX)))
IS_GCC = 1
endif
ifneq (,$(findstring clang++,$(CXX)))
IS_CLANG = 1
endif
ifneq (,$(findstring icpc,$(CXX)))
IS_ICPC = 1
endif

ifdef NOSTRIP
STRIP = echo skip strip
else
STRIP = strip
endif

ifndef QMAKE
QMAKE = qmake
endif

#--------------------------------------------------------
# Define compilation flags to se for each enable feature.
#--------------------------------------------------------

# Flags set to describe a 'prerelease' version.
PRERELEASE = $(shell date +%y%m%d)
PRERELEASE_CFLAGS = -Dgmic_prerelease="\\\"$(PRERELEASE)\\\""

# Minimal set of flags mandatory to compile G'MIC.
MANDATORY_CFLAGS = -Dcimg_date=\\\"\\\" -Dcimg_time=\\\"\\\" -Dcimg_use_zlib -I. $(shell pkg-config --cflags zlib || echo -I$(PREFIX)/$(INCLUDE)) $(PRERELEASE_CFLAGS) $(EXTRA_CFLAGS)
MANDATORY_LIBS = $(shell pkg-config --libs zlib || echo -lz) $(EXTRA_LIBS)

ifndef NO_SRIPDLIB
MANDATORY_CFLAGS += -std=c++11 -pedantic
endif
MANDATORY_LIBS += -L. -L$(PREFIX)/$(LIB) -Wl,-rpath,.
ifdef IS_GCC
MANDATORY_CFLAGS += -Wall -Wextra -Wfatal-errors -Werror=unknown-pragmas -Werror=unused-label
GCC_VER_GTEQ5 = $(shell expr `$(CXX) -dumpversion | cut -f1 -d.` \>= 5)
ifeq ($(GCC_VER_GTEQ5),1)
MANDATORY_CFLAGS += -Wshadow
endif
endif
ifeq ($(OS),Unix)
MANDATORY_CFLAGS += -Dcimg_use_vt100
endif
ifeq ($(OS),Darwin)
ifndef NO_STDLIB
MANDATORY_CFLAGS += -stdlib=libc++
endif
endif
ifeq ($(OS),Windows)
MANDATORY_LIBS += -Wl,--stack,16777216
endif

# Enable optimizations for 'cli'.
FLTO = -flto
OPT_CFLAGS =
ifdef IS_GCC
OPT_CFLAGS += -O3 -mtune=generic # $(FLTO)
OPT_LIBS = # $(FLTO)
endif
ifdef IS_CLANG
OPT_CFLAGS += -O3 -mtune=generic # $(FLTO)
endif
ifdef icpc
OPT_CFLAGS += -fast
endif

# Enable multi-threading support.
PARALLEL_CFLAGS = -Dgmic_is_parallel
ifneq ($(OS),Windows)
PARALLEL_LIBS = -lpthread
endif

# Enable parallelization in CImg, using OpenMP.
# (http://www.openmp.org/)
ifdef IS_CLANG
OPENMP_CFLAGS = -I/usr/lib/gcc/x86_64-linux-gnu/9/include -Dcimg_use_openmp -fopenmp=libomp
OPENMP_LIBS =
else
OPENMP_CFLAGS = -fopenmp -Dcimg_use_openmp
OPENMP_LIBS = -lgomp
endif

# Enable software debugging.
# (Use address sanitizer and thus slows down the code)
ifeq ($(OS),Windows)
DEBUG_CFLAGS = -pedantic -Dcimg_verbosity=3 -Dcimg_strict_warnings -g
else
DEBUG_CFLAGS = -pedantic -Dcimg_verbosity=3 -Dcimg_strict_warnings -g -fsanitize=address # -fsanitize=thread # -fsanitize=undefined # Removed '-ansi' as it (silently) disables c++11.
endif

# Enable extra checking of image buffer pointers.
# (Useful for debug, also slows down the code)
ifneq ($(OS),Darwin)
CHECKIMAGE_CFLAGS = # -Dgmic_check_image
endif

# Enable cancellation points in CImg methods.
ABORT_CFLAGS = -Dcimg_use_abort

# Enable image display, using X11 (Unix).
# (Keep /usr/ dirname here since X11 is located in /usr/ on Mac too)
ifneq ($(OS),Darwin)
X11_CFLAGS = -Dcimg_display=1 $(shell pkg-config --cflags x11 || echo -I/usr/X11R6/include) #-Dcimg_use_xrandr
X11_LIBS = $(shell pkg-config --libs x11 || echo -L/usr/X11R6/lib -lX11) -lpthread # $(shell pkg-config --libs xrandr || echo -lXrandr)
else
ifeq (,$(wildcard /usr/X11))
X11_CFLAGS = -Dcimg_display=0 #-Dcimg_use_xrandr
X11_LIBS = -lpthread # $(shell pkg-config --libs xrandr || echo -lXrandr)
else
X11_CFLAGS = -Dcimg_display=1 $(shell pkg-config --cflags x11 || echo -I/usr/X11/include) #-Dcimg_use_xrandr
X11_LIBS = -L/usr/X11/lib -lX11 -lpthread # $(shell pkg-config --libs xrandr || echo -lXrandr)
endif
endif

# Enable faster X11 display, using XShm extension.
# (ftp://www.x.org/pub/X11R7.7/doc/man/man3/XShm.3.xhtml)
XSHM_CFLAGS = -Dcimg_use_xshm $(shell pkg-config --cflags xcb-shm)
XSHM_LIBS = $(shell pkg-config --libs xcb-shm || echo -L$(PREFIX)/X11R6/lib) -lXext

# Enable image display, using GDI32 (Windows).
GDI32_CFLAGS = -Dcimg_display=2
GDI32_LIBS = -lgdi32

# Enable native support of PNG image files, using the PNG library.
# (http://www.libpng.org/pub/png/libpng.html)
ifneq ($(OS),Darwin)
PNG_CFLAGS = -Dcimg_use_png $(shell pkg-config --cflags libpng)
PNG_LIBS = $(shell pkg-config --libs libpng || echo -lpng -lz)
else
ifeq (,$(wildcard /tmp/skl))
PNG_CFLAGS = -Dcimg_use_png $(shell pkg-config --cflags libpng)
PNG_LIBS = $(shell pkg-config --libs libpng || echo -lpng -lz)
else
PNG_CFLAGS = -Dcimg_use_png
PNG_LIBS = -lpng -lz
endif
endif

# Enable native support of JPEG image files, using the JPEG library.
# (http://libjpeg.sourceforge.net/)
JPEG_CFLAGS = -Dcimg_use_jpeg $(shell pkg-config --cflags libjpeg)
JPEG_LIBS = $(shell pkg-config --libs libjpeg || echo -ljpeg)

# Enable native support of TIFF image files, using the TIFF library.
# (http://www.libtiff.org/)
TIFF_CFLAGS = -Dcimg_use_tiff $(shell pkg-config --cflags libtiff-4)
TIFF_LIBS = $(shell pkg-config --libs libtiff-4 || echo -ltiff)

# Enable native support for loading HEIC/AVIF image files, using the libheif library.
# ( https://github.com/strukturag/libheif )
HEIF_CFLAGS = -Dcimg_use_heif
HEIF_LIBS = -lheif

# Enable native support of MINC2 image files, using the MINC2 library.
# ( http://en.wikibooks.org/wiki/MINC/Reference/MINC2.0_Users_Guide )
MINC2_CFLAGS = -Dcimg_use_minc2 -I${HOME}/local/include
MINC2_LIBS = -lminc_io -lvolume_io2 -lminc2 -lnetcdf -lhdf5 -lz -L${HOME}/local/lib

# Enable native support for downloading files from the network.
# ( http://curl.haxx.se/libcurl/ )
CURL_CFLAGS = -Dcimg_use_curl $(shell pkg-config --cflags libcurl)
CURL_LIBS = $(shell pkg-config --libs libcurl || echo -lcurl)

# Enable native support of webcams and video streaming, using the OpenCV library.
# (https://opencv.org/)
OPENCV_CFLAGS = -Dcimg_use_opencv $(shell pkg-config opencv --cflags) -I/usr/include/opencv -I/usr/include/opencv4
OPENCV_LIBS = $(shell pkg-config opencv --libs || echo -lopencv_core -lopencv_highgui -lopencv_videoio)

# Enable support of most classical image file formats, using the GraphicsMagick++ library.
# (http://www.graphicsmagick.org/Magick++/)
MAGICK_CFLAGS = -Dcimg_use_magick $(shell pkg-config --cflags GraphicsMagick++ || echo -I$(PREFIX)/$(INCLUDE)/GraphicsMagick)
MAGICK_LIBS = $(shell pkg-config --libs GraphicsMagick++ || echo -lGraphicsMagick++)

# Enable native support of EXR image files, using the OpenEXR library.
# (http://www.openexr.com/)
OPENEXR_CFLAGS = -Dcimg_use_openexr $(shell pkg-config --cflags OpenEXR || echo -I$(PREFIX)/$(INCLUDE)/OpenEXR)
OPENEXR_LIBS = $(shell pkg-config --libs OpenEXR || echo -lIlmImf -lHalf)

# Enable Fast Fourier Transforms, using the FFTW3 library.
# (http://www.fftw.org/)
FFTW_CFLAGS = -Dcimg_use_fftw3 $(shell pkg-config --cflags fftw3)
FFTW_LIBS = $(shell pkg-config --libs fftw3 || echo -lfftw3)
FFTW_LIBS += -lfftw3_threads

# Enable native support of the BOARD library.
# (https://github.com/c-koi/libboard)
BOARD_CFLAGS = -Dcimg_use_board
BOARD_LIBS = -lboard

GMIC_STD_CFLAGS = $(MANDATORY_CFLAGS) $(ABORT_CFLAGS) $(PARALLEL_CFLAGS) $(FFTW_CFLAGS) $(CURL_CFLAGS) $(PNG_CFLAGS) $(JPEG_CFLAGS) $(TIFF_CFLAGS)
GMIC_STD_LIBS = $(MANDATORY_LIBS) $(PARALLEL_LIBS) $(FFTW_LIBS) $(CURL_LIBS) $(PNG_LIBS) $(JPEG_LIBS) $(TIFF_LIBS)
ifeq ($(OS),Unix) # Unix
GMIC_STD_CFLAGS += $(OPENMP_CFLAGS) $(X11_CFLAGS) $(OPENEXR_CFLAGS)
GMIC_STD_LIBS += $(OPENMP_LIBS) $(X11_LIBS) $(OPENEXR_LIBS)
else
ifeq ($(OS),Darwin) # MacOSX
GMIC_STD_CFLAGS += $(X11_CFLAGS) $(OPENEXR_CFLAGS)
GMIC_STD_LIBS += $(X11_LIBS) $(OPENEXR_LIBS)
else # Windows
GMIC_STD_CFLAGS += $(OPENMP_CFLAGS) $(GDI32_CFLAGS)
GMIC_STD_LIBS += $(OPENMP_LIBS) $(GDI32_LIBS)
endif
endif

#-------------------------------
# Define main Makefile entries.
#-------------------------------

# Multi-targets.
#---------------
all:
	@echo "**"
ifeq ($(OS),Unix)
	@echo "** Start building G'MIC with default Unix configuration."
else
ifeq ($(OS),Darwin)
	@echo "** Start building G'MIC with default MacOSX configuration."
else
	@echo "** Start building G'MIC with default Windows configuration."
endif
endif
	@echo "**"
	$(MAKE) lib cli gimp gmic_qt libc zart

native:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) -O3 -march=native" "LIBS+=$(GMIC_STD_LIBS)" cli

# libgmic.
#---------
lib:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) $(OPT_CFLAGS)" "LIBS+=$(GMIC_STD_LIBS) $(OPT_LIBS)" _lib

_lib : libgmic.so use_libgmic

libgmic.so: libgmic.o
	ar rcs libgmic.a libgmic.o
ifeq ($(OS),Darwin)
	$(CXX) -shared -std=c++11 -pedantic -o libgmic.so.$(VERSION) libgmic.o $(LIBS)
else
	$(CXX) -shared -std=c++11 -pedantic -Wl,-soname,libgmic.so.$(VERSION1) -o libgmic.so.$(VERSION) libgmic.o $(LIBS)
endif
ifneq ($(OS),Windows)
	ln -fs libgmic.so.$(VERSION) libgmic.so.$(VERSION1)
	ln -fs libgmic.so.$(VERSION1) libgmic.so
endif

libgmic.o: gmic.cpp gmic.h gmic_stdlib.h CImg.h
	$(CXX) -o libgmic.o -c gmic.cpp -Dgmic_core $(PIC) $(CFLAGS)

use_libgmic : use_libgmic.cpp libgmic.so
	$(CXX) -o use_libgmic use_libgmic.cpp -std=c++11 -pedantic -lgmic $(LIBS)

# CLI.
#-----
cli:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) $(OPT_CFLAGS)" "LIBS+=gmic.o $(GMIC_STD_LIBS) $(OPT_LIBS)" gmic.o _cli
	$(STRIP) gmic$(EXE)

dev:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS)" "LIBS+=gmic.o $(GMIC_STD_LIBS)" gmic.o _cli

debug:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) $(DEBUG_CFLAGS)" "LIBS+=gmic.o $(GMIC_STD_LIBS)" gmic.o _cli

valgrind:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) -g" "LIBS+=gmic.o $(GMIC_STD_LIBS)" gmic.o _cli

opencv:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) $(OPT_CFLAGS) $(OPENCV_CFLAGS)" "LIBS+=gmic.o $(GMIC_STD_LIBS) $(OPT_LIBS) $(OPENCV_LIBS)" gmic.o _cli

half:
	$(MAKE) "CFLAGS+=$(MANDATORY_CFLAGS) $(ABORT_CFLAGS) $(PARALLEL_CFLAGS) $(OPENMP_CFLAGS) $(X11_CFLAGS) $(OPENEXR_CFLAGS) -Dcimg_use_half -I/usr/include/OpenEXR -Dgmic_pixel_type=half" \
	        "LIBS+=gmic.o $(MANDATORY_LIBS) $(PARALLEL_LIBS) $(OPENMP_LIBS) $(X11_LIBS) $(OPENEXR_LIBS)" gmic.o _cli

_cli: gmic$(EXE)

gmic$(EXE): gmic_cli.cpp gmic.o
	$(CXX) -o gmic$(EXE) gmic_cli.cpp $(CFLAGS) $(LIBS)

gmic.o: gmic.cpp gmic.h gmic_stdlib.h CImg.h
	$(CXX) -o gmic.o -c gmic.cpp -Dgmic_core $(CFLAGS)

# CLI (shared).
#--------------
cli_shared:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) $(OPT_CFLAGS)" "LIBS+=$(GMIC_STD_LIBS) $(OPT_LIBS) -lgmic" _lib _cli
	$(STRIP) gmic$(EXE)

# CLI (static).
#--------------
GMIC_STATIC_CLI_PATH = $(PREFIX)/$(LIB)/x86_64-linux-gnu
GMIC_STATIC_CLI_EXTRA =
GMIC_STATIC_CLI_CFLAGS = -Dcimg_display=0 $(MANDATORY_CFLAGS) $(OPENMP_CFLAGS) $(ABORT_CFLAGS) $(PARALLEL_CFLAGS) $(PNG_CFLAGS) $(JPEG_CFLAGS) $(FFTW_CFLAGS)
GMIC_STATIC_CLI_LIBS = ./libgmic.o $(PARALLEL_LIBS) \
                       $(GMIC_STATIC_CLI_PATH)/libz.a \
	               $(GMIC_STATIC_CLI_PATH)/libpng.a \
	               $(GMIC_STATIC_CLI_PATH)/libjpeg.a \
	               $(GMIC_STATIC_CLI_PATH)/libfftw3_threads.a $(GMIC_STATIC_CLI_PATH)/libfftw3.a \
	               $(GMIC_STATIC_CLI_EXTRA)

cli_static:
	rm -f libgmic.o
	$(MAKE) "CFLAGS+=$(GMIC_STATIC_CLI_CFLAGS) $(OPT_CFLAGS)" "LIBS+=$(GMIC_STATIC_CLI_LIBS) $(OPT_LIBS)" libgmic.o _cli
	$(STRIP) gmic$(EXE)

# libcgmic.
#----------
libc:
	$(MAKE) "CFLAGS+=$(GMIC_STD_CFLAGS) $(OPT_CFLAGS)" "LIBS+=$(GMIC_STD_LIBS) $(OPT_LIBS)" _lib _libc

_libc: libcgmic_files libcgmic.so use_libcgmic

libcgmic_files:
	@echo "> Retrieve G'MIC libcgmic files..."
	@if [ ! -f gmic_libc.cpp ]; then \
	  if [ -d ../../gmic-community/libcgmic ]; then \
	    ln -fs ../../gmic-community/libcgmic/gmic_libc.cpp .; \
	    ln -fs ../../gmic-community/libcgmic/gmic_libc.h .; \
	    ln -fs ../../gmic-community/libcgmic/use_libcgmic.c .; \
	  else \
            $(WGET) gmic_libc.cpp https://raw.githubusercontent.com/GreycLab/gmic-community/master/libcgmic/gmic_libc.cpp; \
            $(WGET) gmic_libc.h https://raw.githubusercontent.com/GreycLab/gmic-community/master/libcgmic/gmic_libc.h; \
            $(WGET) use_libcgmic.c https://raw.githubusercontent.com/GreycLab/gmic-community/master/libcgmic/use_libcgmic.c; \
            touch gmic_libc.cpp gmic_libc.h use_libcgmic.c; \
	  fi; \
	fi
	@echo " done!"

libcgmic.so: libcgmic.o
	ar rcs libcgmic.a libcgmic.o
ifeq ($(OS),Darwin)
	$(CXX) -shared -o libcgmic.so.$(VERSION) libcgmic.o libgmic.o $(LIBS)
else
ifeq ($(OS),Windows)
	$(CXX) -shared -o libcgmic-$(VERSION1).dll libcgmic.o libgmic.o -Wl,--output-def,libcgmic.def,--out-implib,libcgmic.a $(LIBS)
else
	$(CXX) -shared -Wl,-soname,libcgmic.so.$(VERSION1) -o libcgmic.so.$(VERSION) libcgmic.o libgmic.o $(LIBS)
endif
endif
ifneq ($(OS),Windows)
	ln -fs libcgmic.so.$(VERSION) libcgmic.so.$(VERSION1)
	ln -fs libcgmic.so.$(VERSION1) libcgmic.so
endif

libcgmic.o: gmic_libc.cpp CImg.h gmic_libc.h gmic.h
	$(CXX) -o libcgmic.o -c gmic_libc.cpp -Dgmic_core $(PIC) $(CFLAGS)

use_libcgmic: use_libcgmic.c gmic_libc.h
	$(CC) -std=c99 -o use_libcgmic use_libcgmic.c -lcgmic $(LIBS)

# libcgmic (static).
#-------------------
GMIC_LIBC_STATIC_CFLAGS = $(MANDATORY_CFLAGS) $(ABORT_CFLAGS) $(PARALLEL_CFLAGS) $(FFTW_CFLAGS) $(PNG_CFLAGS) $(JPEG_CFLAGS)
ifeq ($(OS),Window)
GMIC_LIBC_STATIC_LIBS = -Bstatic -lpthread
endif
GMIC_LIBC_STATIC_LIBS += $(MANDATORY_LIBS) $(PARALLEL_LIBS) $(FFTW_LIBS) $(PNG_LIBS) $(JPEG_LIBS)
ifeq ($(OS),Unix) # Unix.
GMIC_LIBC_STATIC_CFLAGS += $(OPENMP_CFLAGS) $(X11_CFLAGS)
GMIC_LIBC_STATIC_LIBS += $(OPENMP_LIBS) $(X11_LIBS)
else
ifeq ($(OS),Darwin) # MacOSX.
GMIC_LIBC_STATIC_CFLAGS += $(X11_CFLAGS)
GMIC_LIBC_STATIC_LIBS += $(X11_LIBS)
else # Windows.
GMIC_LIBC_STATIC_CFLAGS += $(OPENMP_CFLAGS) $(GDI32_CFLAGS)
GMIC_LIBC_STATIC_LIBS += $(OPENMP_LIBS) $(GDI32_LIBS)
endif
endif

libc_static:
ifeq ($(OS),Windows)
	rm libgmic.o
endif
	$(MAKE) "CFLAGS+=$(GMIC_LIBC_STATIC_CFLAGS) $(OPT_CFLAGS)" "LIBS+=$(GMIC_LIBC_STATIC_LIBS) $(OPT_LIBS)" _lib _libc_static

_libc_static: libcgmic_files libcgmic_static.so use_libcgmic_static

libcgmic_static.so: libcgmic_static.o
	ar rcs libcgmic_static.a libcgmic_static.o
ifeq ($(OS),Darwin)
	$(CXX) -shared -static -o libcgmic_static.so.$(VERSION) libcgmic_static.o libgmic.o $(PIC) $(CFLAGS) $(LIBS)
else
ifeq ($(OS),Windows)
	$(CXX) -shared -static -o libcgmic_static-$(VERSION1).dll libcgmic_static.o libgmic.o $(PIC) -Wl,--output-def,libcgmic_static.def,--out-implib,libcgmic_static.a $(CFLAGS) $(LIBS)
else
	$(CXX) -shared -Wl,-soname,libcgmic_static.so.$(VERSION1) -o libcgmic_static.so libcgmic_static.o libgmic.o $(PIC) $(CFLAGS) $(LIBS)
endif
endif
ifneq ($(OS),Windows)
	ln -fs libcgmic_static.so.$(VERSION) libcgmic_static.so.$(VERSION1)
	ln -fs libcgmic_static.so.$(VERSION1) libcgmic_static.so
endif

libcgmic_static.o: gmic_libc.cpp CImg.h gmic_libc.h gmic.h
	$(CXX) -o libcgmic_static.o -c gmic_libc.cpp $(PIC) $(CFLAGS)

use_libcgmic_static: use_libcgmic.c gmic_libc.h
	$(CC) -std=c99 -o use_libcgmic_static use_libcgmic.c -lcgmic_static $(LIBS)

# G'MIC Online.
#--------------
GMICOL_LIB_PATH = /usr/lib/x86_64-linux-gnu/
GMICOL_CFLAGS = $(MANDATORY_CFLAGS) $(PARALLEL_CFLAGS) $(FFTW_CFLAGS) $(PNG_CFLAGS) $(JPEG_CFLAGS) -Dcimg_display=0 # $(OPENMP_CFLAGS)
GMICOL_LIBS = $(MANDATORY_LIBS) $(PARALLEL_LIBS) $(PNG_LIBS) $(JPEG_LIBS) $(GMICOL_LIB_PATH)/libfftw3.a $(GMICOL_LIB_PATH)/libfftw3_threads.a # $(OPENMP_LIBS)

gmicol:
	$(MAKE) "CFLAGS+=$(GMICOL_CFLAGS) $(OPT_CFLAGS)" "LIBS+=$(GMICOL_LIBS) $(OPT_LIBS)" _cli
	$(STRIP) gmic$(EXE)

# G'MIC-Qt.
#----------
gimp: gmic_qt_files CImg.h gmic_stdlib_community.h
	cd ../gmic-qt && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) PRERELEASE=$(PRERELEASE) GMIC_DYNAMIC_LINKING=off HOST=gimp gmic_qt.pro && $(MAKE)
	@echo "Executable 'gmic_gimp_qt' has been successfully compiled in '../gmic-qt/'."

gimp_debug :  gmic_qt_files CImg.h gmic_stdlib_community.h
	cd ../gmic-qt && $(QMAKE) CONFIG+=debug GMIC_PATH=$(SRC_PATH) ASAN=on PRERELEASE=$(PRERELEASE) GMIC_DYNAMIC_LINKING=off HOST=gimp gmic_qt.pro && $(MAKE)
	@echo "Executable 'gmic_gimp_qt' has been successfully compiled in '../gmic-qt/'."

gmic_qt: gmic_qt_files CImg.h gmic_stdlib_community.h
	cd ../gmic-qt && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) PRERELEASE=$(PRERELEASE) GMIC_DYNAMIC_LINKING=off HOST=none gmic_qt.pro && $(MAKE)
	@echo "Executable 'gmic_qt' has been successfully compiled in '../gmic-qt/'."

gmic_qt_debug: gmic_qt_files CImg.h gmic_stdlib_community.h
	cd ../gmic-qt && $(QMAKE) CONFIG+=debug GMIC_PATH=$(SRC_PATH) ASAN=on PRERELEASE=$(PRERELEASE) GMIC_DYNAMIC_LINKING=off HOST=none gmic_qt.pro && $(MAKE)
	@echo "Executable 'gmic_qt' has been successfully compiled in '../gmic-qt/'."

gmic_qt_files:
	@if [ ! -e ../gmic-qt ]; then \
	   curr_pwd=`pwd`; \
	   if [ ! -d ../../gmic-qt ]; then \
	     echo "**"; \
	     echo "** Warning: folder 'gmic-qt' was not found !"; \
	     echo "** It is mandatory to compile the G'MIC plug-in for GIMP ('gmic_gimp_qt'),"; \
	     echo "** as well as the standalone Qt-based interface 'gmic_qt'."; \
             echo "** Trying to retrieve it, with: $$ cd ../../ && git clone https://github.com/c-koi/gmic-qt.git"; \
	     echo "**"; \
	     if cd ../../ && git clone https://github.com/c-koi/gmic-qt.git; then \
	       echo "**"; \
               echo "** Success !"; \
	       echo "**"; \
	     else \
	       echo "**"; \
               echo "** Fatal error: Failed to retrieve folder 'gmic-qt', compilation will probably ends shortly !"; \
	       echo "**"; \
	     fi; \
	   fi; \
	   cd $$curr_pwd; \
	   if [ -d ../../gmic-qt ]; then ln -fs ../gmic-qt ..; fi \
	 fi

# G'MIC-Qt (shared).
#-------------------
gimp_shared: gmic_qt_files CImg.h gmic_stdlib_community.h
	cd ../gmic-qt && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) PRERELEASE=$(PRERELEASE) GMIC_DYNAMIC_LINKING=on HOST=gimp gmic_qt.pro && $(MAKE)
	@echo "Executable 'gmic_gimp_qt' has been successfully compiled in '../gmic-qt/'."

gmic_qt_shared: gmic_qt_files CImg.h gmic_stdlib_community.h
	cd ../gmic-qt && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) PRERELEASE=$(PRERELEASE) GMIC_DYNAMIC_LINKING=on HOST=none gmic_qt.pro && $(MAKE)
	@echo "Executable 'gmic_qt' has been successfully compiled in '../gmic-qt/'."

# ZArt.
#------
zart_shared: zart_files CImg.h gmic_stdlib.h
ifeq ($(OS),Darwin)
	cd ../zart && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) GMIC_DYNAMIC_LINKING=on zart.pro && $(MAKE) && $(STRIP) zart.app/Contents/MacOS/zart
else
	cd ../zart && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) GMIC_DYNAMIC_LINKING=on zart.pro && $(MAKE) && $(STRIP) zart
endif
	@echo "Executable 'zart' has been successfully compiled in '../zart/'."

zart: zart_files CImg.h gmic_stdlib.h
ifeq ($(OS),Darwin)
	cd ../zart && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) zart.pro && $(MAKE) && $(STRIP) zart.app/Contents/MacOS/zart
else
	cd ../zart && $(QMAKE) CONFIG+=release GMIC_PATH=$(SRC_PATH) zart.pro && $(MAKE) && $(STRIP) zart
endif
	@echo "Executable 'zart' has been successfully compiled in '../zart/'."

zart_files:
	@if [ ! -e ../zart ]; then \
	   curr_pwd=`pwd`; \
	   if [ ! -d ../../zart ]; then \
	     echo "**"; \
	     echo "** Warning: folder 'zart' was not found !"; \
	     echo "** It is mandatory to compile 'zart'."; \
             echo "** Trying to retrieve it, with: $$ cd ../../ && git clone https://github.com/c-koi/zart.git"; \
	     echo "**"; \
	     if cd ../../ && git clone https://github.com/c-koi/zart.git; then \
	       echo "**"; \
               echo "** Success !"; \
	       echo "**"; \
	     else \
               echo "**"; \
               echo "** Fatal error: Failed to retrieve folder 'zart', compilation will probably ends shortly !"; \
               echo "**"; \
	     fi; \
	   fi; \
	   cd $$curr_pwd; \
	   if [ -d ../../zart ]; then ln -fs ../zart ..; fi \
	 fi

# Get header files.
#------------------
CImg.h:
	@echo "> Retrieve CImg Library..."
	@if [ -f ../../CImg/CImg.h ]; then \
          if [ ! -f ./CImg.h ]; then ln -s ../../CImg/CImg.h .; fi; \
        elif [ ! -e ./CImg.h ]; then \
          $(WGET) CImg.h https://github.com/GreycLab/CImg/raw/master/CImg.h; \
          touch CImg.h; \
        fi
	@echo " done!"

gmic_stdlib.h: gmic_stdlib.gmic
	@echo "> Retrieve G'MIC Standard Library..."
	@if [ ! -f ./gmic_stdlib.h ]; then \
	  $(WGET) gmic_stdlib.h https://gmic.eu/gmic_stdlib$(VERSION).h; \
	  touch gmic_stdlib.h; \
	fi
	@echo " done!"

do_gmic_stdlib.h:
	@echo "/*\n\
 #\n\
 #  File        : gmic_stdlib.h (v.$(SVERSION))\n\
 #                ( C++ header file )\n\
 #\n\
 #  Description : Raw data arrays encoding the G'MIC standard library\n\
 #                into a compressed form.\n\
 #                This file has been automatically generated by the Makefile entry\n\
 #                'do_gmic_stdlib.h:', from the G'MIC source file 'gmic_stdlib.gmic'.\n\
 #                ( https://gmic.eu )\n\
 #\n\
 #  Copyright   : David Tschumperlé\n\
 #                ( https://tschumperle.users.greyc.fr/ )\n\
 #\n\
 #  Licenses    : This file is 'dual-licensed', you have to choose one\n\
 #      	   of the two licenses below to apply.\n\
 #\n\
 #                CeCILL-C\n\
 #                The CeCILL-C license is close to the GNU LGPL.\n\
 #                ( http://cecill.info/licences/Licence_CeCILL-C_V1-en.html )\n\
 #\n\
 #            or  CeCILL v2.1\n\
 #                The CeCILL license is compatible with the GNU GPL.\n\
 #                ( http://cecill.info/licences/Licence_CeCILL_V2.1-en.html )\n\
 #\n\
 #  This software is governed either by the CeCILL or the CeCILL-C license\n\
 #  under French law and abiding by the rules of distribution of free software.\n\
 #  You can  use, modify and or redistribute the software under the terms of\n\
 #  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA\n\
 #  at the following URL: \"http://cecill.info\".\n\
 #\n\
 #  As a counterpart to the access to the source code and  rights to copy,\n\
 #  modify and redistribute granted by the license, users are provided only\n\
 #  with a limited warranty  and the software's author,  the holder of the\n\
 #  economic rights,  and the successive licensors  have only  limited\n\
 #  liability.\n\
 #\n\
 #  In this respect, the user's attention is drawn to the risks associated\n\
 #  with loading,  using,  modifying and/or developing or reproducing the\n\
 #  software by the user in light of its specific status of free software,\n\
 #  that may mean  that it is complicated to manipulate,  and  that  also\n\
 #  therefore means  that it is reserved for developers  and  experienced\n\
 #  professionals having in-depth computer knowledge. Users are therefore\n\
 #  encouraged to load and test the software's suitability as regards their\n\
 #  requirements in conditions enabling the security of their systems and/or\n\
 #  data to be ensured and,  more generally, to use and operate it in the\n\
 #  same conditions as regards security.\n\
 #\n\
 #  The fact that you are presently reading this means that you have had\n\
 #  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.\n\
 #\n\
*/\n" > gmic_stdlib.h
	@\gmic gmic_stdlib.gmic it gmic_stdlib.gmic parse_gui update "(0)" a y serialize uint8,1,0 o -.h | sed 's/uint8/const unsigned char/' | sed 's/unnamed/gmic/' >> gmic_stdlib.h
	@echo "\nconst unsigned long size_data_gmic = (unsigned long)sizeof(data_gmic);" >> gmic_stdlib.h
	@\gmic v - _upload gmic_stdlib.h
	@cp -f gmic_stdlib.h /tmp/gmic_stdlib$(VERSION).h
	@\gmic v - _upload /tmp/gmic_stdlib$(VERSION).h
	@\gmic v - _upload gmic_stdlib.h

gmic_stdlib_community.h: gmic_stdlib.gmic
	@echo "> Retrieve G'MIC Standard Library + Community Contributions..."
	@if [ ! -f ./gmic_stdlib_community.h ]; then \
	  $(WGET) gmic_stdlib_community.h https://gmic.eu/gmic_stdlib_community$(VERSION).h; \
	  touch gmic_stdlib_community.h; \
	fi
	@echo " done!"

do_gmic_stdlib_community.h:
	@echo "/*\n\
 #\n\
 #  File        : gmic_stdlib_community.h (v.$(SVERSION))\n\
 #                ( C++ header file )\n\
 #\n\
 #  Description : Raw data arrays encoding the G'MIC standard library\n\
 #                plus the community contributions, into a compressed form.\n\
 #                This file has been automatically generated by the Makefile entry\n\
 #                'do_gmic_stdlib_community.h:', from the G'MIC source file 'updateXXX.gmic'.\n\
 #                ( https://gmic.eu )\n\
 #\n\
 #  Copyright   : David Tschumperlé\n\
 #                ( https://tschumperle.users.greyc.fr/ )\n\
 #\n\
 #  Licenses    : This file is 'dual-licensed', you have to choose one\n\
 #      	   of the two licenses below to apply.\n\
 #\n\
 #                CeCILL-C\n\
 #                The CeCILL-C license is close to the GNU LGPL.\n\
 #                ( http://cecill.info/licences/Licence_CeCILL-C_V1-en.html )\n\
 #\n\
 #            or  CeCILL v2.1\n\
 #                The CeCILL license is compatible with the GNU GPL.\n\
 #                ( http://cecill.info/licences/Licence_CeCILL_V2.1-en.html )\n\
 #\n\
 #  This software is governed either by the CeCILL or the CeCILL-C license\n\
 #  under French law and abiding by the rules of distribution of free software.\n\
 #  You can  use, modify and or redistribute the software under the terms of\n\
 #  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA\n\
 #  at the following URL: \"http://cecill.info\".\n\
 #\n\
 #  As a counterpart to the access to the source code and  rights to copy,\n\
 #  modify and redistribute granted by the license, users are provided only\n\
 #  with a limited warranty  and the software's author,  the holder of the\n\
 #  economic rights,  and the successive licensors  have only  limited\n\
 #  liability.\n\
 #\n\
 #  In this respect, the user's attention is drawn to the risks associated\n\
 #  with loading,  using,  modifying and/or developing or reproducing the\n\
 #  software by the user in light of its specific status of free software,\n\
 #  that may mean  that it is complicated to manipulate,  and  that  also\n\
 #  therefore means  that it is reserved for developers  and  experienced\n\
 #  professionals having in-depth computer knowledge. Users are therefore\n\
 #  encouraged to load and test the software's suitability as regards their\n\
 #  requirements in conditions enabling the security of their systems and/or\n\
 #  data to be ensured and,  more generally, to use and operate it in the\n\
 #  same conditions as regards security.\n\
 #\n\
 #  The fact that you are presently reading this means that you have had\n\
 #  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.\n\
 #\n\
*/\n" > gmic_stdlib_community.h
	@\gmic gmic_stdlib.gmic update it '$${_path_rc}'update'$${_version}'.gmic "(0)" a y serialize uint8,1,0 o -.h | sed 's/uint8/const unsigned char/' | sed 's/unnamed/gmic/' >> gmic_stdlib_community.h
	@echo "\nconst unsigned long size_data_gmic = (unsigned long)sizeof(data_gmic);" >> gmic_stdlib_community.h
	@\gmic v - _upload gmic_stdlib_community.h
	@cp -f gmic_stdlib_community.h /tmp/gmic_stdlib_community$(VERSION).h
	@\gmic v - _upload /tmp/gmic_stdlib_community$(VERSION).h
	@\gmic v - _upload gmic_stdlib_community.h

# Bash completion script.
#------------------------
bashcompletion:
	@mkdir -p ../resources
	@\gmic gmic_stdlib.gmic it gmic_stdlib.gmic parse_cli bashcompletion > ../resources/gmic_bashcompletion.sh
	@echo "Bash completion script 'gmic_bashcompletion.sh' has been successfully generated in '../resources/'."

# Man pages.
#-----------
man:
	@mkdir -p ../man
	@\gmic gmic_stdlib.gmic it gmic_stdlib.gmic reference man > ../man/gmic.1
	@gzip -f ../man/gmic.1
	@echo "Man file 'gmic.1.gz' has been successfully generated in '../man/'."

# Install / uninstall / clean.
#-----------------------------
install:
	mkdir -p $(DESTDIR)$(PREFIX)/$(INCLUDE)
	mkdir -p $(DESTDIR)$(PREFIX)/$(BIN)
	mkdir -p $(DESTDIR)$(PLUGINDIR)
	mkdir -p $(DESTDIR)$(PREFIX)/share
	mkdir -p $(DESTDIR)$(PREFIX)/$(LIB)
	@if [ -f gmic ]; then cp -f gmic $(DESTDIR)$(PREFIX)/$(BIN); fi
	@if [ -f gmic_libc.h ]; then cp -f gmic_libc.h $(DESTDIR)$(PREFIX)/$(INCLUDE)/; fi
	@if [ -f ../resources/gmic_cluts.gmz ]; then cp -f ../resources/gmic_cluts.gmz $(DESTDIR)$(PLUGINDIR); fi
	@if [ -f ../resources/gmic_denoise_cnn.gmz ]; then cp -f ../resources/gmic_denoise_cnn.gmz $(DESTDIR)$(PLUGINDIR); fi
	@if [ -f ../zart/zart ]; then cp -f ../zart/zart $(DESTDIR)$(PREFIX)/$(BIN); fi
	@if [ -f ../gmic-qt/gmic_qt ]; then cp -f ../gmic-qt/gmic_qt $(DESTDIR)$(PREFIX)/$(BIN); fi
	@if [ -f ../gmic-qt/gmic_gimp_qt ]; then cp -f ../gmic-qt/gmic_gimp_qt $(DESTDIR)$(PLUGINDIR); fi
	@if [ -f libgmic.so ]; then \
	  cp -f gmic.h $(DESTDIR)$(PREFIX)/$(INCLUDE); \
          cp -f libgmic.so $(DESTDIR)$(PREFIX)/$(LIB)/libgmic.so.$(VERSION); \
	  ln -fs libgmic.so.$(VERSION) $(DESTDIR)$(PREFIX)/$(LIB)/libgmic.so.$(VERSION1); \
	  ln -fs libgmic.so.$(VERSION1) $(DESTDIR)$(PREFIX)/$(LIB)/libgmic.so; \
	fi
	@if [ -f libcgmic.so ]; then \
	  cp -f gmic_libc.h $(DESTDIR)$(PREFIX)/$(INCLUDE); \
	  cp -f libcgmic.so $(DESTDIR)$(PREFIX)/$(LIB)/libcgmic.so.$(VERSION); \
          ln -fs libcgmic.so.$(VERSION) $(DESTDIR)$(PREFIX)/$(LIB)/libcgmic.so.$(VERSION1); \
          ln -fs libcgmic.so.$(VERSION1) $(DESTDIR)$(PREFIX)/$(LIB)/libcgmic.so; \
        fi
	@if [ -f ../man/gmic.1.gz ]; then \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/man/man1; \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/man/fr/man1; \
	  cp -f ../man/gmic.1.gz $(DESTDIR)$(PREFIX)/share/man/man1/gmic.1.gz; \
	  cp -f ../man/gmic.1.gz $(DESTDIR)$(PREFIX)/share/man/fr/man1/gmic.1.gz; \
	fi
	@if [ -f ../resources/gmic_bashcompletion.sh ]; then \
	  mkdir -p $(DESTDIR)$(BASH_COMPLETIONSDIR); \
          cp -f ../resources/gmic_bashcompletion.sh $(DESTDIR)$(BASH_COMPLETIONSDIR)/gmic; \
        fi
	@if [ -f ../zart/zart ]; then \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/applications; \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/icons/hicolor/48x48/apps; \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/icons/hicolor/scalable/apps; \
	  cp -f ../zart/zart.desktop $(DESTDIR)$(PREFIX)/share/applications; \
	  cp -f ../zart/icons/48-zart.png $(DESTDIR)$(PREFIX)/share/icons/hicolor/48x48/apps/zart.png; \
	  cp -f ../zart/icons/zart.svg $(DESTDIR)$(PREFIX)/share/icons/hicolor/scalable/apps; \
	fi
	@if [ -f ../gmic-qt/gmic_qt ]; then \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/applications; \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/icons/hicolor/48x48/apps; \
	  mkdir -p $(DESTDIR)$(PREFIX)/share/icons/hicolor/scalable/apps; \
	  cp -f ../gmic-qt/gmic_qt.desktop $(DESTDIR)$(PREFIX)/share/applications; \
	  cp -f ../gmic-qt/icons/application/48-gmic_qt.png $(DESTDIR)$(PREFIX)/share/icons/hicolor/48x48/apps/gmic_qt.png; \
	  cp -f ../gmic-qt/icons/application/gmic_qt.svg $(DESTDIR)$(PREFIX)/share/icons/hicolor/scalable/apps; \
	fi

uninstall:
	rm -f $(DESTDIR)$(PLUGINDIR)/gmic_gimp_qt
	rm -f $(DESTDIR)$(PREFIX)/$(BIN)/gmic
	rm -f $(DESTDIR)$(PREFIX)/$(BIN)/zart
	rm -f $(DESTDIR)$(PREFIX)/$(BIN)/gmic_qt
	rm -f $(DESTDIR)$(PREFIX)/$(INCLUDE)/gmic.h
	rm -f $(DESTDIR)$(PREFIX)/$(LIB)/libgmic.so.$(VERSION)
	rm -f $(DESTDIR)$(PREFIX)/$(LIB)/libgmic.so.$(VERSION1)
	rm -f $(DESTDIR)$(PREFIX)/$(LIB)/libgmic.so
	rm -f $(DESTDIR)$(PREFIX)/$(LIB)/libcgmic.so.$(VERSION)
	rm -f $(DESTDIR)$(PREFIX)/$(LIB)/libcgmic.so.$(VERSION1)
	rm -f $(DESTDIR)$(PREFIX)/$(LIB)/libcgmic.so
	rm -rf $(DESTDIR)$(PREFIX)/share/doc/gmic/
	rm -f $(DESTDIR)$(PREFIX)/share/man/man1/gmic.1.gz
	rm -f $(DESTDIR)$(PREFIX)/share/man/fr/man1/gmic.1.gz
	rm -f $(DESTDIR)$(PREFIX)/share/applications/zart.desktop
	rm -f $(DESTDIR)$(PREFIX)/share/icons/hicolor/48x48/apps/zart.png
	rm -f $(DESTDIR)$(PREFIX)/share/icons/hicolor/scalable/apps/zart.svg
	rm -f $(DESTDIR)$(PREFIX)/share/applications/gmic_qt.desktop
	rm -f $(DESTDIR)$(PREFIX)/share/icons/hicolor/48x48/apps/gmic_qt.png
	rm -f $(DESTDIR)$(PREFIX)/share/icons/hicolor/scalable/apps/gmic_qt.svg

distclean: clean

clean:
	rm -rf CImg.h gmic_stdlib.h gmic_stdlib_community.h gmic$(EXE) use_libgmic$(EXE) use_libcgmic$(EXE) use_libcgmic_static$(EXE) gmic*.o libgmic* libcgmic* *~
	@if [ -f ../zart/Makefile ]; then cd ../zart && $(MAKE) clean; fi
	@if [ -h ../zart ]; then rm -f ../zart; fi
	@if [ -f ../gmic-qt/Makefile ]; then cd ../gmic-qt && $(MAKE) clean; fi
	@if [ -h ../gmic-qt ]; then rm -f ../gmic-qt; fi
	@if [ -h gmic_libc.h ]; then rm -f gmic_libc.h; fi
	@if [ -h gmic_libc.cpp ]; then rm -f gmic_libc.cpp; fi
	@if [ -h use_libcgmic.c ]; then rm -f use_libcgmic.c; fi

# End of Makefile.
