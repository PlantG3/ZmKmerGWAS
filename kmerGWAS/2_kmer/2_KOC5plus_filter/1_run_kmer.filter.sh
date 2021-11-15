### merge:
outdir=KOC5plus
mkdir $outdir
for txt in ../../1-kmers/counts/*txt; do
        echo $txt
        geno=$(echo $txt | sed 's/.*\///g' | sed 's/.txt//g')
        out=$outdir/"SB.".$geno.sbatch
        echo "#!/bin/bash -l" > $out
        echo "#SBATCH --mem-per-cpu=4G" >> $out
        echo "#SBATCH --time=1-00:00:00" >> $out
        echo "#SBATCH --cpus-per-task=1" >> $out
        echo "awk '\$2>=5' "$txt" > "$outdir"/"$geno".txt" >> $out
        sbatch $out
done