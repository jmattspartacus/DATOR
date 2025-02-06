LIBDIR=$(shell pwd)/lib/
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

libReader :
	cd src/Reader && $(MAKE)

libGRETINA : libReader
	cd src/GRETINA && $(MAKE)

libORRUBA : libReader
	cd src/ORRUBA && $(MAKE)

libS800 : libReader
	cd src/S800 && $(MAKE)

LDFMerge : LDFMerge.c
	$(C) -std=c99 -O3 -o LDFMerge LDFMerge.c -lz

LDFConvert : LDFConvert.c
	$(C) -std=c99 -O3 -o LDFConvert LDFConvert.c -lz

install : libReader libGRETINA libORRUBA libS800
	cd src/Reader && $(MAKE) install
	cd src/GRETINA && $(MAKE) install
	cd src/ORRUBA && $(MAKE) install
	cd src/S800 && $(MAKE) install
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/bin
	cp LDFMerge $(INSTALLDIR)/bin
	cp LDFConvert $(INSTALLDIR)/bin	

clean:
	rm lib/*.so
	rm LDFMerge LDFConvert
