#!/bin/bash
dataPath=<path to trimmed reads>
perl kmer.sbatch.pl \
	--mem 3G --time 0-12:00:00 \
	--indir $dataPath --outdir counts \
	--fq1feature .R1.pair.fq.gz --fq2feature .R2.pair.fq.gz \
	--kmersize 31 \
	--mincount 1 \
	--hashsize "500M" \
	--fa2txt_script "./fasta2txt.pl" \
	--jellyfish "/path to jellyfish/jellyfish" \
	--checkscript \
	--threads 8

