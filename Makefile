
BUILDDIR=$(shell pwd)/build
LIBDIR=$(BUILDDIR)/lib/
BINDIR=$(BUILDDIR)/bin/
BUILDINCLUDEDIR=$(BUILDDIR)/include/DATOR

INSTALLDIR=$(HOME)/.local

CC = g++
CFLAGS = `root-config --cflags` -O3 -g -fPIC
LIBS = -Wl,--no-as-needed -lz
ROOTLIBS = `root-config --libs --glibs` -Wl,--no-as-needed -lMathMore
C = gcc

export CC
export CFLAGS
export LIBS
export LIBDIR
export INSTALLDIR

all: libReader libGRETINA libORRUBA libS800 LDFMerge LDFConvert


prepdir :
	mkdir -p $(BINDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BUILDINCLUDEDIR)
	cd src && find . -type f -name "*.hh" -exec cp --parents {} ../build/include/DATOR/ \;

libReader : prepdir
	cd src/Reader && $(MAKE)

libGRETINA : libReader
	cd src/GRETINA && $(MAKE)

libORRUBA : libReader
	cd src/ORRUBA && $(MAKE)

libS800 : libReader
	cd src/S800 && $(MAKE)


LDFMerge : prepdir LDFMerge.c
	$(C) -std=c99 -O3 -o $(BINDIR)LDFMerge LDFMerge.c -lz

LDFConvert : prepdir LDFConvert.c
	$(C) -std=c99 -O3 -o $(BINDIR)LDFConvert LDFConvert.c -lz

all : libReader libGRETINA libORRUBA libS800 LDFMerge LDFConvert

install : all
	cp -r $(BUILDDIR)/* $(INSTALLDIR)

clean:
	rm -rf build
