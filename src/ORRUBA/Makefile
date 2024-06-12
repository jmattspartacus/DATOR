all: libORRUBA

libORRUBA : ORRUBA.cc ORRUBA.hh
	$(CC) $(CFLAGS) -I../ -o $(LIBDIR)/libORRUBA.so -shared ORRUBA.cc $(LIBS) -L../../lib/ -lReader

install : $(LIBDIR)libORRUBA.so
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/include/DATOR/ORRUBA/	
	cp $(LIBDIR)/libORRUBA.so $(INSTALLDIR)/lib/
	cp *.hh $(INSTALLDIR)/include/DATOR/ORRUBA/
