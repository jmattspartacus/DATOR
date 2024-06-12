all: libReader

libReader : Reader.cc Reader.hh RawData.hh
	$(CC) $(CFLAGS) --shared -o $(LIBDIR)/libReader.so Reader.cc $(LIBS)

install:
	mkdir -p $(INSTALLDIR)/lib/
	mkdir -p $(INSTALLDIR)/include/DATOR/Reader/
	cp $(LIBDIR)/libReader.so $(INSTALLDIR)/lib/
	cp *.hh $(INSTALLDIR)/include/DATOR/Reader/

