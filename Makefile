
PROGRAMS = mrfast_output_to_mipcounts detail_mip_targets #call_mip_srgap2_cn
BINARIES = $(addprefix bin/,$(PROGRAMS))

CC = gcc
DEBUG = -g
CFLAGS = -Wall -lz

all : $(BINARIES)

.PHONY : clean

bin/% : src/%.c
	$(CC) $(CFLAGS) $^ -o $@

clean :
	rm -f $(BINARIES)
