#!/bin/bash
#SBATCH --array=1-518
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=./slurmLog/slurm_%A_%a_%j.out

id=$SLURM_ARRAY_TASK_ID  

##########################################################################
### input: subject to changes
##########################################################################
# R script
rkgwas=1m-kgwas_copynumber.R
# koc directory
kmdir=<path of k-mer counts directory>
# koc mean depth
mean_depth=4.1
# trait
trait=<path of trait file (.csv format)>
# PCA (pop structure)
pc=<path of population structure file>
# max pvalue
maxpval=1e-6
# output directory
outdir=`echo $trait | sed 's/.*\///g'`
##########################################################################

##########################################################################
# run
##########################################################################
if [ ! -d $outdir ]; then
        mkdir $outdir
fi

# run R:
module load R
Rscript $rkgwas $id $mean_depth $kmdir $trait $pc $outdir $maxpval

