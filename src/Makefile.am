# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libimperial.a libimperial.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libimperial.so ${prefix}/lib
	cp libimperial.a ${prefix}/lib
	cp imperial.h ${prefix}/include

libimperial.a: imperial_static.o 
	$(AR) rcs $@ $^

libimperial.so: imperial.o
	$(FC) -shared $(AM_FCFLAGS) -o libimperial.so $^ $(AM_LDFLAGS)

imperial.o: imperial.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
imperial_static.o: imperial.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

