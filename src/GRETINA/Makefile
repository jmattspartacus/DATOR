all: libGRETINA

libGRETINA : Gretina.cc Gretina.hh GretinaConf.cc GretinaConf.hh GretinaHit.cc GretinaHit.hh Gamma.hh Gamma.cc
	$(CC) $(CFLAGS) -I../ -o $(LIBDIR)libGRETINA.so -shared Gretina.cc GretinaConf.cc GretinaHit.cc Gamma.cc $(LIBS) -L$(LIBDIR) -lReader

install: $(LIBDIR)/libGRETINA.so
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/include/DATOR/GRETINA/	
	cp $(LIBDIR)/libGRETINA.so $(INSTALLDIR)/lib
	cp *.hh $(INSTALLDIR)/include/DATOR/GRETINA/
