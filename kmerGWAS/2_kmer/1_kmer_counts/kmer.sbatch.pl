use warnings;
use Getopt::Long;

my ($fqfilecol, $mem, $time, $jellyfish, $indir, $db, $outdir, $help);
my ($kmersize, $mincount, $hashsize, $fa2txt_script);
my ($fq1feature, $fq2feature, $threads, $checkscript);
my $sleep_sec = 0;
my $result = &GetOptions("mem=s" => \$mem,
            "time=s" => \$time,
            "jellyfish=s" => \$jellyfish,
			"indir=s" => \$indir,
			"outdir=s" => \$outdir,
			"fq1feature=s" => \$fq1feature,
			"fq2feature=s" => \$fq2feature,
			"threads=i" => \$threads,
			"kmersize=i" => \$kmersize,
			"mincount=i" => \$mincount,
			"hashsize=s" => \$hashsize,
			"fa2txt_script=s" => \$fa2txt_script,
			"sleep=i" => \$sleep_sec,
			"checkscript" => \$checkscript,
			"help|h" => \$help
);

# print help information if errors occur:
if ($help) {
	&errINF;
	exit;
}

$mem = "36G" if (!defined $mem);
$time = "48:00:00" if (!defined $time);
$fq1feature = ".R1.pair.fq" if (!defined $fq1feature);
$fq2feature = ".R2.pair.fq" if (!defined $fq2feature);
$outdir = "." if (!defined $outdir);
if (!-d $outdir) {
	`mkdir $outdir`;
}
$threads = 8 if (!defined $threads);
$kmersize = 25 if (!defined $kmersize);
$mincount = 1 if (!defined $mincount);
$hashsize = "500M" if (!defined $hashsize);
$fa2txt_script = "/homes/liu3zhen/scripts/kmer/fasta2txt.pl" if (!defined $fa2txt_script);
$jellyfish = "/homes/liu3zhen/local/bin/jellyfish" if (!defined $jellyfish);

open (IN,"ls \"$indir\" -1 |");
while (<IN>) {
	chomp;
	my $fqfile = $_;
	if ($fqfile =~ $fq1feature) {
		my $sample = $fqfile;
		$sample =~ s/.*\///g;
		$sample =~ s/$fq1feature//g;
		$sample = $outdir."\/".$sample;
		my $fq1 = $fqfile;
		my $fq2 = $fq1;
		$fq2 =~ s/$fq1feature/$fq2feature/g;
		print "$fqfile\n";
		print "$sample\n";
		my $outfile = $sample.".sbatch";
		my $out = $sample.".sam";
		
		my $infq1 = $indir."/".$fq1;
		if ($fq1 =~ /gz$/) {
			$infq1 = "<(gunzip -c ".$indir."/".$fq1.")";
		}	
		
		my $infq2 = $indir."/".$fq2; 
		if ($fq2 =~ /gz$/) {
			$infq2 = "<(gunzip -c ".$indir."/".$fq2.")";
		}
		
		# sleep to avoid the problem for job assignments
		`sleep $sleep_sec`;
		
		open(OUT, ">", $outfile) || die;
		print OUT "#!/bin/bash -l\n";
		print OUT "#SBATCH --mem-per-cpu=$mem\n";
		print OUT "#SBATCH --time=$time\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks-per-node=$threads\n";
		print OUT "$jellyfish count -m $kmersize -t $threads -o $sample.jf -L $mincount -s $hashsize -C $infq1 $infq2\n";
		print OUT "$jellyfish dump $sample.jf > $sample.fa\n";
		print OUT "$fa2txt_script --kcutoff $mincount $sample\.fa\n";
		print OUT "rm $sample.jf; rm $sample.fa\n";
		close OUT;
		my $sbatch_cmd = sprintf("sbatch %s", $outfile);
		#print "$sbatch_cmd\n";
		if ($checkscript) {
			print STDERR "NO running; SBATCH scripts were generated\n";
		} else {
			system($sbatch_cmd);
		}
	}
}

close IN;

sub errINF {
	# to be added (SL)
}
