cap3=/homes/liu3zhen/software/cap3/CAP3/cap3
pval_cutoff=2.118e-08
koc_inpath=/bulk/liu3zhen/research/projects/kmerGWAS/main2/1a_oil/1a1_kgCn/1a1_kgCn
koc_path=`realpath $koc_inpath`
prefix=oil.cn.kgwas
log=${prefix}.log
min_mapping_score=30
final_prefix=1a2o_oil.cn.asm_B97

B97ref=/homes/liu3zhen/references/NAM1.0/1-db/bwa/B97-1.0.fasta

# create output directory
if [ ! -d $final_prefix ]; then
	mkdir $final_prefix
fi

cd $final_prefix

### 1. extract k-mers
echo "extract k-mers"
### positive correlated k-mers:
awk -v p=$pval_cutoff '$2>0 && $3<p' ${koc_path}/*KOC.txt > ${prefix}.positive.kmers.KOC
awk '{ print ">k"NR"_"$3"\n"$1 }' ${prefix}.positive.kmers.KOC > ${prefix}.positive.kmers.KOC.fasta

### negative correlated k-mers:
awk -v p=$pval_cutoff '$2<0 && $3<p' ${koc_path}/*KOC.txt > ${prefix}.negative.kmers.KOC
awk '{ print ">k"NR"_"$3"\n"$1 }' ${prefix}.negative.kmers.KOC > ${prefix}.negative.kmers.KOC.fasta


### 2. assembly with cap3
echo "assemble k-mers"
for kmfas in *.kmers.KOC.fasta; do
	$cap3 $kmfas -i 25 -j 35 -s 300 -o 20 -p 94 >$log >&1
	#perl ~/scripts/cap3/ace.ctg.reads.pl ${kmfas}.cap.ace 1>/dev/null 2>${kmfas}.cap.contig.pval.means
	#perl ~/scripts/cap3/fasta.name.replace.pl ${kmfas}.cap.contigs ${kmfas}.cap.contig.pval.means > ${kmfas}.cap.contigs.rename
	perl /homes/liu3zhen/scripts/cap3/ace.ctg.reads.minP.pl ${kmfas}.cap.ace 1>/dev/null 2>${kmfas}.cap.contig.min.pval
	perl /homes/liu3zhen/scripts/cap3/fasta.name.replace.pl ${kmfas}.cap.contigs ${kmfas}.cap.contig.min.pval > ${kmfas}.cap.contigs.rename
	cat ${kmfas}.cap.contigs.rename ${kmfas}.cap.singlets > ${kmfas}.cap.fasta
	perl /homes/liu3zhen/corescripts/fasta/fasta.size.filter.pl --min 0 --max 69 ${kmfas}.cap.fasta > ${kmfas}.cap.lt70.fasta
	perl /homes/liu3zhen/corescripts/fasta/fasta.size.filter.pl --min 70 ${kmfas}.cap.fasta > ${kmfas}.cap.gte70.fasta
done

### 3. align to reference
echo "align seq"
# small sequences:
for smfas in *lt70.fasta; do
	/homes/liu3zhen/software/bwa/bwa-0.7.17/bwa aln $B97ref $smfas > ${smfas}.B97.sai
	/homes/liu3zhen/software/bwa/bwa-0.7.17/bwa samse -n 5 -f ${smfas}.B97.sam $B97ref ${smfas}.B97.sai $smfas
	awk -v s=$min_mapping_score '$5 > s' ${smfas}.B97.sam > ${smfas}.filter.B97.sam

done

#KC.kmergwas.positive.kmers.KOC.fasta.cap.lt70.fasta.B73.sai


# large sequences:
for lgfas in *gte70.fasta; do
	/homes/liu3zhen/software/bwa/bwa-0.7.17/bwa mem $B97ref $lgfas > ${lgfas}.B97.sam
	awk -v s=$min_mapping_score '$5 > s' ${lgfas}.B97.sam > ${lgfas}.filter.B97.sam

done

### 4. merge data
echo "merge sam"
cat *negative*filter.B97.sam | grep "^@" -v | sort -k3n,3 -k4n,4 > ${final_prefix}.negative.mapping.B97.sam
cat *positive*filter.B97.sam | grep "^@" -v | sort -k3n,3 -k4n,4 > ${final_prefix}.positive.mapping.B97.sam

cd ..

# copy significant k-mers
head -n 1 ${koc_path}/*_1.normKOC*kgwas.KOC.txt > tmp.header 
cat tmp.header $final_prefix/${prefix}.positive.kmers.KOC > ${final_prefix}.positive.kmers.KOC
cat tmp.header $final_prefix/${prefix}.negative.kmers.KOC > ${final_prefix}.negative.kmers.KOC
cat tmp.header $final_prefix/${prefix}.positive.kmers.KOC $final_prefix/${prefix}.negative.kmers.KOC > ${final_prefix}.associated.kmers.KOC


