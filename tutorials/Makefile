INSTALLDIR=$(HOME)/.local/

CC = g++
CFLAGS = `root-config --cflags` -O3 -g -fPIC
LIBS = -Wl,--no-as-needed -lz
ROOTLIBS = `root-config --libs --glibs` -Wl,--no-as-needed -lMathMore

all : BasicSort

BasicSort : BasicSort.cc libReader libGRETINA libORRUBA
	$(CC) $(CFLAGS) -I$(INSTALLDIR)/include/DATOR/ -o BasicDatorSort BasicSort.cc $(ROOTLIBS) $(LIBS) -L$(LIBDIR) -lReader -lGRETINA -lORRUBA
