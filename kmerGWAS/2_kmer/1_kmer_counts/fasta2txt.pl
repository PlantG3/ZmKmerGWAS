#!/usr/bin/perl -w
# count fa to table
# Sanzhen Liu
# June 18th 2014

use strict;
use warnings;
use Getopt::Long;

my ($kmer_cutoff, $help);
my $result = &GetOptions(
	"kcutoff=i" => \$kmer_cutoff,
	"help" => \$help);

if ($help) {
	&errINF;
}

my $count;
my $seq;

foreach my $eachfile (@ARGV) {
	my $sample = $eachfile;
	$sample =~ s/\.fa$|\.fas$|\.fasta$//g;
	my $outfile = $sample.".txt";
	open(OUT, ">$outfile") || die;
	print OUT "Kmer\t$sample\n";

	open(IN, $eachfile) || die;
	while (<IN>) {
		chomp;
		if (/^>(.+)/) {
			### output
			if (defined $count and $count >= $kmer_cutoff) {
				print OUT "$seq\t$count\n";
			}
			$count = $1;
			$seq = '';
		} else {
			$seq .= $_;
		}
	}
	### output last element
	if ($count >= $kmer_cutoff) {
		print OUT "$seq\t$count\n";
	}
	close IN;
	close OUT;
}

