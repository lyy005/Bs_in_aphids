#!/usr/bin/perl
use strict;

die "perl $0 [collinearity file] [gff file] [spp1 name as in gff] [spp2 name as in gff] [output]\n" unless @ARGV == 5;

my $spp1 = $ARGV[2];
my $spp2 = $ARGV[3];

my %SMI_gff_s;
my %API_gff_s;
my %SMI_gff_e;
my %API_gff_e;

# save gene coordinates in gff
open (GFF, "$ARGV[1]") or die "$ARGV[1] $!\n";
while(<GFF>){
	chomp;
	my $line = $_;
	my @line = split(/\t/,$line);
	if($line[0] =~ /$spp1/){
		#print "$line[0]\n";
		$SMI_gff_s{$line[1]} = $line[2];
		$SMI_gff_e{$line[1]} = $line[3];
	}
	if($line[0] =~ /$spp2/){
		#print "$line\n"; 
		$API_gff_s{$line[1]} = $line[2];
		$API_gff_e{$line[1]} = $line[3];
	}
}
close GFF;

# recode the start and end of each colinear block
open (CO, "$ARGV[0]") or die "$ARGV[0] $!\n";
open(OUT,">$ARGV[4]") or die "$ARGV[4] $!\n";
open(LOG,">$ARGV[4].log") or die "$ARGV[4].log $!\n";

my $ID;
my %SMI;
my %API;

while(<CO>){
	chomp;
	my $line = $_;

	if($line =~ /^##/ and $line !~ /Alignment/){
		next;
	}
	if($line =~ /Alignment/){
		$ID = $line;
		print LOG "$ID\n";
	}elsif($ID){
		my @line = split(/\t/,$line);

		print LOG "$line[1] from $SMI_gff_s{$line[1]} to $SMI_gff_e{$line[1]} +++ $line[2] from $API_gff_s{$line[2]} to $API_gff_e{$line[2]}\n";

		push @{$SMI{$ID}},$line[1]; # SMI as spp1, also comes first (column 2) in the collinear block
		push @{$API{$ID}},$line[2]; # API as spp2, also comes second (column 3) in the collinear block
	}
}
close CO;

#print @{$SMI{'## Alignment 51: score=6355.0 e_value=0 N=129 SMI2&API4 minus'}},"\n";

my %final;
foreach my $each (sort keys %API){
	my $API_s = ${$API{$each}}[0];
	my $API_e = ${$API{$each}}[-1];
	my $SMI_s = ${$SMI{$each}}[0];
	my $SMI_e = ${$SMI{$each}}[-1];

	my $SMI_block_s = $SMI_gff_s{$SMI_s};
	my $SMI_block_e = $SMI_gff_e{$SMI_e};
	my $API_block_s = $API_gff_s{$API_s};
	my $API_block_e = $API_gff_e{$API_e};

	my @arr = split(/\s/,$each);
	my @aaa = split(/&/,$arr[6]);
	#$final{"$aaa[0]\t$API_block_s\t$API_block_e\t$aaa[1]\t$SMI_block_s\t$SMI_block_e"} = $each;
	print OUT "$aaa[0]\t$SMI_block_s\t$SMI_block_e\t$aaa[1]\t$API_block_s\t$API_block_e\n";
}

#foreach my $each (sort(keys %final)){
	#print OUT "$final{$each}\n$each\n";
#	print OUT "$each\n";
#}
close OUT;

print "Outputs:
 $ARGV[4]
 $ARGV[4].log\n";
