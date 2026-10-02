#!/usr/bin/perl -w
use strict;

die "perl $0 [mitoscaf.fa.gbf.gene.fasta] [output dir]\n" unless @ARGV == 2;

my %hash;
open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";

$/=">";
<IN>;
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\n+/;
	my $id = shift @line;
	my $seq = join "", @line;
	my @id = split /\;/, $id;

	if(defined($id[1])){
		if($hash{$id[1]}){
			$hash{$id[1]} .= ">$id[0]\n$seq\n";
		}else{
			$hash{$id[1]}  = ">$id[0]\n$seq\n";
		}
	}else{
		print "Format error: $id\n";
	}
}
close IN;
$/="\n";

print "Output: \n";
foreach my $k (sort keys %hash){
	print " $ARGV[1]/$k.nt.fasta\n";
	open (OUT, ">$ARGV[1]/$k.nt.fasta") or die "$ARGV[1]/$k.nt.fasta $!\n";
	print OUT "$hash{$k}";
	close OUT;
}
