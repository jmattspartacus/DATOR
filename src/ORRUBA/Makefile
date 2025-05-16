all: libORRUBA

libORRUBA : ORRUBA.cc ORRUBA.hh ORRUBAAuxiliary.cc ORRUBAAuxiliary.hh SX3.cc SX3.hh BB10.cc BB10.hh QQQ5.cc QQQ5.hh SingleParticle.hh TDC.hh TDC.cc
	$(CC) $(CFLAGS) -I../ -o $(LIBDIR)/libORRUBA.so -shared ORRUBA.cc QQQ5.cc SX3.cc BB10.cc ORRUBAAuxiliary.cc TDC.cc $(LIBS) -L$(LIBDIR) -lReader
