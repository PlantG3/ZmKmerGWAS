#!/usr/bin/perl -w
# ===============================================================
# kmer.split.pl
# usage: perl kmer.split.pl [kmer files]
# ===============================================================

use strict;
use warnings;
use File::Temp;

my @ntc;
my @nt = ("A", "T", "C", "G");
foreach my $nt1 (@nt) {
	foreach my $nt2 (@nt) {
		foreach my $nt3 (@nt) {
			my $ntc = $nt1.$nt2.$nt3;
			push(@ntc, $ntc);
			system(sprintf("%s%s", "mkdir ",$ntc));
		}
	}
}

### for each Kmer file
foreach my $kmerfile (@ARGV) {
	for (my $i = 0; $i <= $#ntc; $i++) {
		my $prop = $i/($#ntc + 1);
		print STDERR "$prop has been done for $kmerfile.\n";
		my $comb = $ntc[$i];
		my $match = "^".$comb;
		my $kmerfilename = $kmerfile;
		$kmerfilename =~ s/.*\///g;
		my $out = "./".$comb."/".$kmerfilename.".$comb";
		#printf("%s%s %s%s%s%s", "grep -P -e ", $comb, $kmerfile, " > ", $out);
		system(sprintf("%s%s %s%s%s", "grep -P -e ", $match, $kmerfile, " > ", $out));
	}
}
