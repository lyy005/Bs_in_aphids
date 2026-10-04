#!/usr/bin/perl -w
use strict;

die "perl $0 [GTF file] [order file] [output file]

* version 0.9 - sort scaffolds based on an order file
 # sort numerically based on scaffold IDs
 # only tested on Braker Tsebra gtf files\n" unless @ARGV == 3;

my %hash;
my @order;
open (OR, "$ARGV[1]") or die "$ARGV[1] $!\n";
while(<OR>){
	chomp;
	my @line = split;
	push @order, $line[0];
}
close OR;

my $order = join " ", @order;
print "Scaffold order is: 
$order\n"; 

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
while(<IN>){
	my $all = $_;
	my @line = split /\s+/;
	my $id = $line[0];
	
	if($hash{$id}){
		$hash{$id} .= $all;
	}else{
		$hash{$id} = $all;
	}
}
close IN;

# sort scaffolds by the numerical parts of their IDs
# scaffold_1
#foreach my $key (sort {($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0]} keys %hash){
foreach my $key (@order){
	if($hash{$key}){
		print OUT "$hash{$key}";
	}else{
		print "No genes on $key\n";
	}
}
close OUT;
