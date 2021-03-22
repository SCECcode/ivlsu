# Autoconf/automake file

bin_PROGRAMS = test_imperial

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_LDFLAGS = ${LDFLAGS}

objects = test.o
TARGETS = $(bin_PROGRAMS)

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/tests
	cp test_imperial ${prefix}/tests

test_imperial: $(objects)
	$(CC) -o $@ $^ $(AM_CFLAGS) -L../src -limperial $(AM_LDFLAGS)

$(objects): %.o: %.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS) -I../src/ $(AM_CFLAGS)
