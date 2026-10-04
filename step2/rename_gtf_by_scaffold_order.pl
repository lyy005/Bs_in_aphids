#!/usr/bin/perl -w
use strict;

die "perl $0 [GTF file] [output file]\n # rename gtf gene names based on the numerical order of scaffold IDs\n # only tested on Braker Tsebra gtf files\n" unless @ARGV == 2;

my %hash;
my $prefix = "g";
my $count = 0;
my $new_gene_id = $prefix.$count;
my $old_gene_id;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OT1, ">$ARGV[1]") or die "$ARGV[1] $!\n";
open (OT2, ">$ARGV[1].cast") or die "$ARGV[1].cast $!\n";
while(<IN>){
	my $line = $_;
	my @line = split /\s+/;

	if($line[2] eq "gene"){
		$count ++;
		$new_gene_id = $prefix.$count;

		$old_gene_id = $line[-1];
		print OT2 "$old_gene_id\t$new_gene_id\t$line";
	}

	$line =~ s/$old_gene_id/$new_gene_id/g;
	print OT1 $line;
}
close IN;

print "Output files:\n $ARGV[1]\n $ARGV[1].cast\n";
