all: libReader

libReader : Reader.cc Reader.hh RawData.hh
	$(CC) $(CFLAGS) --shared -o $(LIBDIR)/libReader.so Reader.cc $(LIBS)


