#!/usr/bin/perl -w
use strict;

die "perl $0 [genome.fa.telo] [length file] [output file]\n
 * v0.92 - updated the format for Table S1
 * v0.91 - check if more than 2 telomeres are found
 * v0.90 - updated the script so it reports scaffolds with telomeres on both ends" unless @ARGV == 3;

my %hash;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\s+/;

	my $region = $line[1]."\t".$line[2];

	# count as left telomere if telomere starts from 0 bp in the genome
	if($line[1] == 0){
		# check if the left telomere already exists
		if($hash{$line[0]}{"left"}){
			print "Left telomere already exists: $all\n";
		}else{
			$hash{$line[0]}{"left"} = $region;
		}
	# else count as rigth telomere
	}else{
		# check if the right telomere already exists
		if($hash{$line[0]}{"right"}){
			print "Right telomere already exists: $all\n";
		}else{
			$hash{$line[0]}{"right"} = $region;
		}
	}
}
close IN;

open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
open (LS, "$ARGV[1]") or die "$ARGV[1] $!\n";
print OUT "scaf\ttotal_len\tleft_s\tleft_e\tright_s\tright_e\n";
while(<LS>){
	my @line = split;
	print OUT "$line[0]\t$line[1]\t";
	if($hash{$line[0]}{"left"}){
		print OUT "$hash{$line[0]}{'left'}\t";
	}else{
		print OUT "NA\tNA\t";
	}

	if($hash{$line[0]}{"right"}){
		print OUT "$hash{$line[0]}{'right'}\n";
	}else{
		print OUT "NA\tNA\n";
	}
}
print "Output: $ARGV[2]\n";
