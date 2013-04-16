
PROGRAMS = mrfast_output_to_mipcounts #detail_mip_targets call_mip_srgap2_cn
BINARIES = $(addprefix bin/,$(PROGRAMS))

CC = gcc
DEBUG = -g
CFLAGS = -Wall -lz -I./iniparser/src
LFLAGS = -L./iniparser -liniparser

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
