#!/usr/bin/perl -w
use strict;

die "perl $0 [combined.collinearity.bed.orient.sv_classified.tsv] [sp1.chr.list] [sp2.chr.list] [output file]
Updates: 
 v20260411 - Some aphids have multiple Bs. Therefore, they have more chrs but not affecting syntenies. So this version, we only use chrs with synteny to calculate the inter-chr rates
 v20260407 - Use the spp with fewer chrs as the ref. The order of the chr list has to be the same order as the sp1 and sp2 cols of the input; In the output, the spp with fewer chrs is in the front
 v20251031 - Combined both outputs into one files
 v20251027 - Added the function to produce both by length and by number of the rearrangements\n" unless @ARGV == 4;

my %chr_pair;
my %sp1_chr;

my %chr_pair_num;
my %sp1_chr_num;

my ($sp1, $sp2);
my $reverse = 0; # whether to reverse the order of the spp in the table
open (SP1, "$ARGV[1]") or die "$ARGV[1] $!\n";
while(<SP1>){
	chomp;
	$sp1 ++;
}
close SP1;

open (SP2, "$ARGV[2]") or die "$ARGV[2] $!\n";
while(<SP2>){
	chomp;
	$sp2 ++;
}
close SP2;

# count the number of chrs with synteny
#print "Chr pairs with synteny: \n";
#print "sp1 - sp2\n";

my (%sp1, %sp2);
#my ($sp1_name, $sp2_name);
open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OT2, ">$ARGV[3].log") or die "$ARGV[3].log $!\n";

my $header = <IN>;
while(<IN>){
	chomp;
	my $line = $_;
	my @line = split /\s+/, $line;

	#print "$line[0] - $line[3]\n";

	$sp1{$line[0]} ++;
	$sp2{$line[3]} ++;

#	$sp1_name = $line[0];
#	$sp2_name = $line[3];
}
close IN;

my ($sp1_syn, $sp2_syn);
my (@sp1_syn, @sp2_syn);
$sp1_syn = keys %sp1;
$sp2_syn = keys %sp2;

@sp1_syn = sort keys %sp1;
@sp2_syn = sort keys %sp2;

print "sp1: @sp1_syn\nsp2: @sp2_syn\n";

# Ref sp must have fewer chrs. By default, sp1 is the ref. If sp1 actually has more chrs, it needs to be reversed.
if($sp1_syn <= $sp2_syn){
	$reverse = 0;
	print "*** Not reversing ***\n";
	print "Total chrs:\n Sp1: $ARGV[1] Number of chrs: $sp1\n Sp2: $ARGV[2] Number of chrs: $sp2\n";
	print "Chrs with synteny:\n Sp1: $ARGV[1] Number of chrs: $sp1_syn\n Sp2: $ARGV[2] Number of chrs: $sp2_syn\n";

	print OT2 "*** Not reversing ***\n";
	print OT2 "Total chrs:\n Sp1: $ARGV[1] Number of chrs: $sp1\n Sp2: $ARGV[2] Number of chrs: $sp2\n";
	print OT2 "Chrs with synteny:\n Sp1: $ARGV[1] Number of chrs: $sp1_syn\n Sp2: $ARGV[2] Number of chrs: $sp2_syn\n";
}else{
	$reverse = 1;
	print "** Reversing the sp1 sp2 order as sp2 has fewer chrs **\n";
	print "Total chrs:\n Sp1: $ARGV[1] Number of chrs: $sp1\n Sp2: $ARGV[2] Number of chrs: $sp2\n";
	print "Chrs with synteny:\n Sp1: $ARGV[1] Number of chrs: $sp1_syn\n Sp2: $ARGV[2] Number of chrs: $sp2_syn\n";

	print OT2 "** Reversing the sp1 sp2 order as sp2 has fewer chrs **\n";
	print OT2 "Total chrs:\n Sp1: $ARGV[1] Number of chrs: $sp1\n Sp2: $ARGV[2] Number of chrs: $sp2\n";
	print OT2 "Chrs with synteny:\n Sp1: $ARGV[1] Number of chrs: $sp1_syn\n Sp2: $ARGV[2] Number of chrs: $sp2_syn\n";
}

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OT, ">$ARGV[3]") or die "$ARGV[3] $!\n";
#open (OT2, ">$ARGV[1].len") or die "$ARGV[1].len $!\n";

$header = <IN>;
chomp $header;
#print OUT "$header\tinter_chr\n";

while(<IN>){
	chomp;
	my $line = $_;
	my @line = split /\s+/, $line;

	if($reverse){
		$sp1_chr{$line[3]} += $line[7];
		$sp1_chr_num{$line[3]} ++;

		my $pair = $line[3]."_".$line[0];
		$chr_pair{$pair} += $line[7];
		$chr_pair_num{$pair} ++;
		
	}else{
		$sp1_chr{$line[0]} += $line[6];
		$sp1_chr_num{$line[0]} ++;

		my $pair = $line[0]."_".$line[3];
		$chr_pair{$pair} += $line[6];
		$chr_pair_num{$pair} ++;
	}
}
close IN;

print OT "chr_pairs\tchr_pair_syn_num\tchr_pair_syn_len\tref_chr\tref_num\tnum_percent\tref_len\tlen_percent\n";

foreach my $k (sort keys %chr_pair){
	my @k = split/\_/, $k;

	my $len_perc = $chr_pair{$k}/$sp1_chr{$k[0]};
	my $num_perc = $chr_pair_num{$k}/$sp1_chr_num{$k[0]};

	print OT "$k\t$chr_pair_num{$k}\t$chr_pair{$k}\t$k[0]\t$sp1_chr_num{$k[0]}\t$num_perc\t$sp1_chr{$k[0]}\t$len_perc\n";
}
print "Output: 
 Sp1: $ARGV[1]
 Sp2: $ARGV[2] 
 Summary of inter-chr rearrangements: $ARGV[3]
 Headers are: 
   col_1 - chr_pairs: chr pairs that share syntenic blocks
   col_2 - chr_pair_syn_num: number of syntenic blocks between the chr pair
   col_3 - chr_pair_syn_len: total length (bp) of syntenic blocks between the chr pair
   col_4 - ref_chr: ref chr name (sp1 chr name)
   col_5 - ref_num: total number of syntenic blocks on that ref chr (chr of sp1)
   col_6 - num_percent: percentage of the number of sytenic blocks in the reference chr, i.e. col_2 / col_4
   col_7 - ref_len: total length of syntenic blocks on that ref chr (chr of sp1)
   col_8 - len_percent: percentage of the length of this chr pair compared to the ref chr, i.e. col_3 / col_5
";
