#!/usr/bin/perl -w
use strict;

die "perl $0 [pep fasta with stop codon in asterisks] [output file]\n" unless @ARGV == 2;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[1]") or die "$ARGV[1] $!\n";

$/=">";
<IN>;
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\n+/;
	my $id = shift @line;
	my $seq = join "", @line;

	$seq=~s/\*$//; # remove the last stop codon

	if($seq=~/\*/){
		print "In frame stop codon found: $id\n";
	}else{

		if($seq=~/\-/){
			print "Incomplete codon found: $id\n";
			$seq=~s/\-//g;
		}
		print OUT ">$id\n$seq\n";
	}
}
close IN;
print "Output: $ARGV[1]\n";
