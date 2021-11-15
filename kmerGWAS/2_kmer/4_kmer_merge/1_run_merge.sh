### merge:
for acgt in /bulk/liu3zhen/public_databases/maize/maize282/8-kmers.merge/1-split/[ACGT]*; do
	echo $acgt
	l3=$(echo $acgt | sed 's/.*\///g')
	out="SB."$l3.sbatch
	echo "#!/bin/bash -l" > $out
	echo "#SBATCH --mem-per-cpu=124G" >> $out
	echo "#SBATCH --time=0-23:00:00" >> $out
	echo "#SBATCH --cpus-per-task=1" >> $out
	echo "#SBATCH --partition=ksu-plantpath-liu3zhen.q,batch.q,killable.q" >> $out
	echo "cat "$acgt"/* | awk 'seen[\$1]++==10' > "$l3".merge" >> $out
	sbatch $out
done

