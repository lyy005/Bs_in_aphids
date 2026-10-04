#!/usr/bin/perl -w
use strict;

die "perl $0 [CDS fasta] [output file]\n # sort based on gene IDs\n" unless @ARGV == 2;

my %hash;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
$/=">";
<IN>;
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\n+/;
	my $id = shift @line;
	my $seq = join "", @line;
#	my @names = split /\_/, $line[0];
#	if($hash{$names}){
	my @id = split /\s+/, $id;

	if($hash{$id[0]}){
		$/="\n";
#		print OUT ">$id\n$seq\n";
		print OUT ">$all";
		$/=">";
	}

}
close IN;
