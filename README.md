MIPs/SUNK genotyper
===================

Dependencies
------------

  * gcc 4.7.0
  * gsl 1.15
  * iniparser 3.1 (http://ndevilla.free.fr/iniparser/iniparser-3.1.tar.gz)

Installation and usage
----------------------

See documentation in docs/mips_sunk_genotyping.html

Genotyping example
------------------

[Download input files and scripts](http://eichlerlab.gs.washington.edu/mips_cnv_typer/mips_cnv_typer_example.tar.gz) for a working example of the genotyper program.

On a Linux-style command line, you can start the example with the following commands.

    tar zxvf mips_cnv_typer_example.tar.gz
    cd mips_cnv_typer_example/
    . config.sh
    export GENOTYPER_DIR=/path/to/mips_cnv_typer/bin
    make

This example requires [mrFAST](http://mrfast.sourceforge.net/) to be in your PATH and expects the compiled mips_cnv_typer binaries to be present in `$GENOTYPER_DIR`.
