# Autoconf/automake file

bin_PROGRAMS = test_ivlsu run_test_ivlsu.sh

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${PROJ_INCLUDES} -I../src
AM_LDFLAGS = ${LDFLAGS} ${PROJ_LDFLAGS} -L../src -livlsu -lm

objects = test_api.o
TARGETS = $(bin_PROGRAMS)

.PHONY = run_test

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/test
	cp test_ivlsu ${prefix}/test

test_ivlsu: $(objects)
	$(CC) -o $@ $^ $(AM_LDFLAGS) 

run_unit : test_ivlsu
	./test_ivlsu

$(objects): %.o: %.c
	$(CC) -c $@ $^ $(AM_CFLAGS) 

clean :
	rm -rf *~ *.o test_ivlsu

