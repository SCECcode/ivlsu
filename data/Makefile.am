# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
	if [ ! -d "ivlsu" ]; then \
		./make_data_files.py; \
	fi

install:
	mkdir -p ${prefix}/data
	cp ./config ${prefix}/data
	cp -R ./ivlsu ${prefix}/data
	chmod og+r ${prefix}/data/config
	chmod -R og+r ${prefix}/data/ivlsu
	chmod og+x ${prefix}/data/ivlsu
