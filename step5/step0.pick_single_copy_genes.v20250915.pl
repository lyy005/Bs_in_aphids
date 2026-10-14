#!/usr/bin/perl -w
use strict;

die "perl $0 [Orthogroups.txt] [sp1] [sp2] [output file]
 * Added a function so can manually input the species names (IDs)\n" unless @ARGV == 4;

my %hash;
open (TAB, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[3]") or die "$ARGV[3] $!\n";

my $sp1 = $ARGV[1];
my $sp2 = $ARGV[2];

# skip header
my $header = <TAB>;
print "Header line: \n$header";

while(<TAB>){
	chomp;
	my $line = $_;
	my @line = split /\s+/, $line; 

	my $c_sp1 = ($line =~ s/$sp1\_/$sp1\_/g);  # gene IDs always start with spp\_geneID
	my $c_sp2 = ($line =~ s/$sp2\_/$sp2\_/g);

	my ($g_sp1, $g_sp2);
	if(($c_sp1 == 1)&&($c_sp2 == 1)){
		if($line =~ /($sp1\_\S+)/){
			$g_sp1 = $1;
		}
	
		if($line =~ /($sp2\_\S+)/){
			$g_sp2 = $1;
		}
		#print OUT "$line[0]\t$g_acpi\t$g_hoco\n";
		print OUT "$g_sp1\t$g_sp2\n";
	}
}
close TAB;

print "Output: $ARGV[3]\n";
