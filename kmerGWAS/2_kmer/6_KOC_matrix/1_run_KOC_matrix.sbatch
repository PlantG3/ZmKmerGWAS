#!/bin/bash
#SBATCH --array=1-624
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=20-00:00:00

id=$SLURM_ARRAY_TASK_ID
kmer=../4-partition/maize282.k31.min10m5_$id
counts=../../1-kmers/counts/*txt
#../1-KOC5plus/KOC5plus/*txt
outDir=KOC31
if [ ! -d $outDir ]; then
        mkdir $outDir
fi
out=maize282.k31.min10m5_${id}.counts
log=maize282.k31.min10m5_${id}.log
perl select.kmer.merge.batch.pl $kmer $counts > $outDir/$out