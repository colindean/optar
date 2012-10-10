LDFLAGS=-lm
CFLAGS=-O3 -Wall -Wuninitialized -fomit-frame-pointer -funroll-loops \
	-fstrength-reduce -DNODEBUG `libpng-config --I_opts`

all: optar unoptar

install:
	install optar /usr/local/bin/
	install unoptar /usr/local/bin
	install pgm2ps /usr/local/bin

uninstall:
	rm /usr/local/bin/optar
	rm /usr/local/bin/unoptar
	rm /usr/local/bin/pgm2ps

clean:
	rm -f optar unoptar golay golay_codes.c *.o

common.o: common.c optar.h
	gcc -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

parity.o: parity.c
	gcc -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

optar.o: optar.c optar.h font.h parity.h
	gcc -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

golay_codes.o: golay_codes.c
	gcc -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

golay.o: golay.c parity.h
	gcc -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

unoptar.o: unoptar.c optar.h parity.h
	gcc -c -I/usr/local/include/libpng $(CPPFLAGS) $(CFLAGS) -o $@ $<

optar: optar.o common.o golay_codes.o parity.o
	gcc $(LDFLAGS) -o $@ $^

golay_codes.c: golay
	./$< > $@

golay: golay.o parity.o
	gcc $(LDFLAGS) -o $@ $^

unoptar: unoptar.o common.o golay_codes.o parity.o
	gcc $(LDFLAGS) -o $@ -L/usr/local/lib -lpng -lz $^
