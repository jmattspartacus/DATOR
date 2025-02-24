all: libORRUBA

libORRUBA : ORRUBA.cc ORRUBA.hh ORRUBAAuxiliary.cc ORRUBAAuxiliary.hh SX3.cc SX3.hh BB10.cc BB10.hh QQQ5.cc QQQ5.hh SingleParticle.hh
	$(CC) $(CFLAGS) -I../ -o $(LIBDIR)/libORRUBA.so -shared ORRUBA.cc QQQ5.cc SX3.cc BB10.cc ORRUBAAuxiliary.cc $(LIBS) -L../../lib/ -lReader

install : $(LIBDIR)libORRUBA.so
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/include/DATOR/ORRUBA/	
	cp $(LIBDIR)/libORRUBA.so $(INSTALLDIR)/lib/
	cp *.hh $(INSTALLDIR)/include/DATOR/ORRUBA/
