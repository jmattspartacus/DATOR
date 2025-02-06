all: libORRUBA

libORRUBA : ORRUBA.cc ORRUBA.hh ORRUBAAuxiliary.cc ORRUBAAuxiliary.hh
	$(CC) $(CFLAGS) -I../ -o $(LIBDIR)/libORRUBA.so -shared ORRUBA.cc ORRUBAAuxiliary.cc $(LIBS) -L../../lib/ -lReader

install : $(LIBDIR)libORRUBA.so
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/include/DATOR/ORRUBA/	
	cp $(LIBDIR)/libORRUBA.so $(INSTALLDIR)/lib/
	cp *.hh $(INSTALLDIR)/include/DATOR/ORRUBA/
