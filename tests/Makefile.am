# Autoconf/automake file

bin_PROGRAMS = test_ivlsu run_test_ivlsu.sh

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${PROJ_INCLUDES} -I../src
AM_LDFLAGS = ${LDFLAGS} ${PROJ_LDFLAGS} -L../src -livlsu -lm

objects = test.o
TARGETS = $(bin_PROGRAMS)

.PHONY = run_test

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/tests
	cp test_ivlsu ${prefix}/tests
	cp run_test_ivlsu.sh ${prefix}/tests

test_ivlsu$(EXEEXT): $(objects)
	$(CC) -g -o $@ $^ $(AM_CFLAGS) $(AM_LDFLAGS) 

$(objects): %.o: %.c
	$(CC) -g -o $@ -c $^ $(AM_CFLAGS) -I${prefix}/include

run_test_ivlsu.sh$(EXEEXT): 

clean :
	rm -rf *~ *.o test_ivlsu

run_test : run_test_ivlsu.sh
	./run_test_ivlsu.sh

