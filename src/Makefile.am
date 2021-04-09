# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libivlsu.a libivlsu.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libivlsu.so ${prefix}/lib
	cp libivlsu.a ${prefix}/lib
	cp ivlsu.h ${prefix}/include

libivlsu.a: ivlsu_static.o 
	$(AR) rcs $@ $^

libivlsu.so: ivlsu.o
	$(CC) -shared $(AM_CFLAGS) -o libivlsu.so $^ $(AM_LDFLAGS)

ivlsu.o: ivlsu.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
ivlsu_static.o: ivlsu.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

