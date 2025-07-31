all: libS800

libS800 : S800.cc S800.hh
	$(CC) $(CFLAGS) -I../ -o $(LIBDIR)/libS800.so -shared S800.cc S800Confs.cc S800Physics.cc $(LIBS) -L$(LIBDIR)/ -lReader

install : $(LIBDIR)libS800.so
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/include/DATOR/S800/	
	cp $(LIBDIR)/libS800.so $(INSTALLDIR)/lib/
	cp *.hh $(INSTALLDIR)/include/DATOR/S800/
