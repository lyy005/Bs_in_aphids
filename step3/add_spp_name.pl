#!usr/bin/perl -w
use strict;

die "Usage: perl $0 [.fa] [prefix of reads] [out.fa]\n" unless (@ARGV == 3);
open (FA, $ARGV[0]) or die "$ARGV[0] $!\n";
open OUT, ">$ARGV[2]" or die "$ARGV[2] $!\n";
open CAST, ">$ARGV[2].cast" or die "$ARGV[2] $!\n";
my $pre = $ARGV[1];

$/=">";
my $null = <FA>;
while(<FA>){
	chomp;
	my @line = split /\n/;
	my $name = shift @line;
	my $seq = join "", @line;
	print OUT ">$pre\_$name\n$seq\n";
	print CAST "$name\t$pre\_$name\n";
}

close FA;
close OUT;
close CAST;
print "Output: 
 $ARGV[2]
 $ARGV[2].cast\n";
