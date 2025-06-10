LIBDIR=$(shell pwd)/build/lib/
BINDIR=$(shell pwd)/build/bin/
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

prepdirs : 
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)

libReader : prepdirs
	cd src/Reader && $(MAKE)

libGRETINA : libReader
	cd src/GRETINA && $(MAKE)

libORRUBA : libReader
	cd src/ORRUBA && $(MAKE)

libS800 : libReader
	cd src/S800 && $(MAKE)

LDFMerge : prepdirs LDFMerge.c
	$(C) -std=c99 -O3 -g -o $(BINDIR)/LDFMerge LDFMerge.c -lz

LDFConvert : prepdirs LDFConvert.c
	$(C) -std=c99 -O3 -g -o $(BINDIR)/LDFConvert LDFConvert.c -lz

install : libReader libGRETINA libORRUBA libS800
	cp src/Reader/*.hh $(INSTALLDIR)/include/DATOR/Reader
	cp src/GRETINA/*.hh $(INSTALLDIR)/include/DATOR/GRETINA
	cp src/S800/*.hh $(INSTALLDIR)/include/DATOR/S800
	cp src/ORRUBA/*.hh $(INSTALLDIR)/include/DATOR/ORRUBA
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/bin
	cp $(BINDIR)* $(INSTALLDIR)/bin/
	cp $(LIBDIR)* $(INSTALLDIR)/lib/

clean:
	rm -rf build
