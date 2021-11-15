### merge:
for txt in ../0-counts/*txt; do
	echo $txt
	geno=$(echo $txt | sed 's/.*\///g' | sed 's/.txt//g')
	out="SB."$geno.sbatch
	echo "#!/bin/bash -l" > $out
	echo "#SBATCH --mem-per-cpu=4G" >> $out
	echo "#SBATCH --time=0-23:00:00" >> $out
	echo "#SBATCH --cpus-per-task=4" >> $out
	echo "perl kmer.split.pl "$txt >> $out
	sbatch $out
done

