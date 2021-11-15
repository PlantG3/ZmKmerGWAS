#!/usr/bin/perl -w
# ===============================================================
# select.kmer.merge.pl
# usage: perl select.kmer.merge.pl <Input1> <Input2>
# Input1: kmer file 1: ONLY one col with Kmer
# Input2: kmer file 2; ONLY 2 cols with Kmer at the 1st col
# ===============================================================
use strict;
use warnings;
use File::Temp;

my (%kmer_hash01, %kmer_hash02);
my (%genocounts);
my $kmer01 = $ARGV[0];
my @kmer02 = @ARGV[1..$#ARGV];

### kmer counts:
my $total_rows = 0;
my $row = 0;
open(IN, $kmer01) || die;
while (<IN>) {
	chomp;
	$row++;
	my @kmercount = split(/\t/, $_);
	my $km = $kmercount[0];
	$genocounts{$row} = $km;
	$kmer_hash01{$km}++;
	$total_rows++;
}
close IN;

### files of lookup kmers
my (@counts, @genolist);
foreach my $kmer02 (@kmer02) {
	open(IN, $kmer02) || die;
	while (<IN>) {
		chomp;
		my ($km, $kmc) = split(/\t/, $_); # kmer and counts
		if (exists $kmer_hash01{$km}) {
			$kmer_hash02{$km} = $kmc;
		}
	}
	
	### output:
	my $row = 0;
	my $geno = $kmer02;
	$geno =~ s/.*\///g;
	$geno =~ s/\.txt//g;
	push(@genolist, $geno);
	open(IN, $kmer01) || die;
	while (<IN>) {
		$row++;
		chomp;
		my @kma = split(/\t/, $_);
		my $km = $kma[0];
		my $current_count = 0;
		if (exists $kmer_hash02{$km}) {
			$current_count = $kmer_hash02{$km};
		}
		my $accum_counts = $genocounts{$row};
		$accum_counts .= "\t".$current_count;
		$genocounts{$row} = $accum_counts;
	}
	print STDERR "finished $geno\n";
	### empty count hash
	%kmer_hash02 = ();
}

### output
# header
print "Kmer";
foreach my $genoname (@genolist) {
	print "\t$genoname";
}
print "\n";

# counts:
for (my $i = 1; $i <= $total_rows; $i++) {
	print "$genocounts{$i}\n";
}

###################################################
