#!/bin/bash
#SBATCH --array=1-624
#SBATCH --output=slurms/slurm_%A_%a_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=1-00:00:00

id=$SLURM_ARRAY_TASK_ID
koc=/path to raw KOC matrix/maize282.k31.min10m5_${id}.counts
outDir=normKOC31

if [ ! -d $outDir ]; then
	mkdir $outDir
fi

module load R
Rscript 1i-rm.lines.and.norm.R $koc $outDir 

