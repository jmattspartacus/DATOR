all: libS800

libS800 : S800.cc S800.hh
	$(CC) $(CFLAGS) -I../ -o $(LIBDIR)/libS800.so -shared S800.cc S800Confs.cc S800Physics.cc $(LIBS) -L$(LIBDIR) -lReader

