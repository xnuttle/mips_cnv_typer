
PROGRAMS = detail_mip_targets mrfast_output_to_mipcounts call_mip_cn
BINARIES = $(addprefix bin/,$(PROGRAMS))

MOD_GSL_DIR = /net/gs/vol3/software/modules-sw/gsl/1.15/Linux/RHEL6/x86_64
CC = gcc
DEBUG = -g
CFLAGS = -Wall -lz -I./iniparser/src -I$(MOD_GSL_DIR)/include
LFLAGS = -L./iniparser -liniparser -L$(MOD_GSL_DIR)/lib -lgsl -lgslcblas

all : libraries $(BINARIES)

.PHONY : clean libraries

bin/% : src/%.c
	mkdir -p bin
	$(CC) $^ $(CFLAGS) $(LFLAGS) -o $@

libraries : iniparser
	cd iniparser && $(MAKE)

iniparser : iniparser-3.1.tar.gz
	tar zxvf $<

iniparser-3.1.tar.gz :
	wget http://ndevilla.free.fr/iniparser/iniparser-3.1.tar.gz

clean :
	rm -f $(BINARIES)
	cd lib/iniparser && $(MAKE) veryclean
