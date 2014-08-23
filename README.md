
========================================================================================================================

IRBIS: a systematic search for conserved complementarity

This repository contains a pipeline for conserved complementarity that comes along with
the publication in RNA 2014 (PMID:25142064)

A short manual is included (see Manual/ folder) to give minimal description of the IRBIS pipeline.
Currently the pipeline is designed to replicate the results of the publication PMID:25142064, but
it can be re-configured to perform user-specific RNA structure seqrches. Questions shall be addressed 
to Dmitri Pervouchine (dp at crg dot eu).


========================================================================================================================

INSTALLATION

Check out this repository from github. The pipeline is based on GNU make that is used for two 
purposes: (i) compiling the code and (ii) data processing and workflow

Enter the folder and run sequentially

	make all
	make -f database.mk all
	make -f metacalc.mk all

The first step ('make all') creates binaries and configures makefiles for the following steps. 

The second step ('make -f database.mk all') does automatic data download from public repositories, including
genomic sequences, pairwise sequence alignments, and annotations. The download location is specified
in config.dat file. Note that $FIXEDPATH is the path to the database directory. You might want to skip 
the download step (at least for genomes) if you already have the corresponding files.

The third step ('make -f metacalc.mk all') creates metadata in the $METADATA directory and runs IRBIS pipeline.
Output files will be located on $OUTDIR.


IRBIS depends on g++, perl, latex + tikz + pgfplots (pgfplots is optional), and muscle (if muscle software is
not in the system, you can obtain it by 'make muscle'). Important: check the $muscle variable in Perl/align.pm

========================================================================================================================

For more detailed information on examples, benchmarks, and running the pipeline in custom modes
please refer to the manual.
