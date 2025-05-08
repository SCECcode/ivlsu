# Autoconf/automake file

bin_PROGRAMS = test_ivlsu run_test_ivlsu.sh

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} -I../src
AM_LDFLAGS = ${LDFLAGS} -L../src -livlsu

objects = test.o
TARGETS = $(bin_PROGRAMS)

.PHONY = run_test

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/test
	cp test_ivlsu ${prefix}/test
	cp run_test_ivlsu.sh ${prefix}/test

test_ivlsu$(EXEEXT): $(objects)
	$(CC) -g -o $@ $^ $(AM_CFLAGS) $(AM_LDFLAGS) 

$(objects): %.o: %.c
	$(CC) -g -o $@ -c $^ $(AM_CFLAGS) -I${prefix}/include

run_test_ivlsu.sh$(EXEEXT): 

clean :
	rm -rf *~ *.o test_ivlsu

run_test : run_test_ivlsu.sh
	./run_test_ivlsu.sh

