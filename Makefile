
PROGRAMS = mrfast_output_to_mipcounts #call_mip_srgap2_cn detail_mip_targets
BINARIES = $(addprefix bin/,$(PROGRAMS))

CC = gcc
DEBUG = -g
CFLAGS = -Wall -lz

all : $(BINARIES)

bin/% : src/%.c
	$(CC) $(CFLAGS) $^ -o $@
