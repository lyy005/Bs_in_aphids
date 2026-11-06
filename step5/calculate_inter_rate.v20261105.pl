#!/usr/bin/perl -w
use strict;

die "perl $0 [inter_sv.list] [pairwise_time.tsv] [output file]
 v20261105 - add inter-rate per MY = the number of rearrangements / Mb of syntenic blocks / MY
 Rearrangement rates are calculated as the number of inter-chr rearrangements / total alignment (syntenic block) length of sp1 / divergence time\n" unless @ARGV == 3;

open (IN, "$ARGV[1]") or die "$ARGV[1] $!\n";
open (OT2, ">$ARGV[2].log") or die "$ARGV[2].log $!\n";

my $header = <IN>;
print $header;

my %time;
while(<IN>){
	chomp;
	my $line = $_;
	my @line = split/\s+/;

	my @spp = ($line[0], $line[1]);
	my $spp = join "_", @spp;

	if($time{$spp}){
		if($time{$spp} == $line[3]){

		}else{
			print "Time is not consistant in the table: $spp $time{$spp} $line[3]\n";
		}
	}else{
		print OT2 "Divergence time: $spp $line[3]\n";
		$time{$spp} = $line[3];
	}
}
close IN;

open (LS, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";

#print OUT "spp\ttime\tinter_chr_counts\trates\tauto_auto_counts\tauto_auto_rates\tX_X_counts\tX_X_rates\tauto_X_counts\tauto_X_rates\n";
print OUT "spp\ttime\tInterchr_block_rate\tInterchr_block_rate_AtoA\tInterchr_block_rate_XtoX\tInterchr_block_rate_AtoX\tAvg_number_of_breaks_perMb\tAvg_number_of_breaks_perMb_AtoA\tAvg_number_of_breaks_perMb_XtoX\tAvg_number_of_breaks_perMb_AtoX\tinter_rate\tInterchr_block_rate_perMY\tInterchr_block_rate_AtoA_perMY\tInterchr_block_rate_XtoX_perMY\tInterchr_block_rate_AtoX_perMY\tAvg_number_of_breaks_perMb_perMY\tAvg_number_of_breaks_perMb_AtoA_perMY\tAvg_number_of_breaks_perMb_XtoX_perMY\tAvg_number_of_breaks_perMb_AtoX_perMY\tinter_rate_perMY\n";

while(<LS>){
	chomp;
	my $line = $_;
	my @line = split/\//, $line;
	my @spp = split/\_/, $line[0];
	print "####
# Species pairs: $spp[0] - $spp[1]
####\n";

# Interchr_block_rate	Interchr_block_rate_AtoA	Interchr_block_rate_XtoX	Interchr_block_rate_AtoX	Avg_number_of_breaks_perMb	Avg_number_of_breaks_perMb_AtoA	Avg_number_of_breaks_perMb_XtoX	Avg_number_of_breaks_perMb_AtoX
# 0	0	0	0	0	0	0	0
## Interchr_block_rate = off-diagnal (bp)/ Total syntenic blocks (bp)
## Avg_number_of_breaks_perMb = Off-diagnal breaks / Off-diagnal length

	open (SV, "$line") or die "$line $!\n";
	my $header = <SV>;
	print $header;

	my $rates;
	my @rates;
	my $inter_rate;
	my $inter_rate_AtoA;
	my $inter_rate_XtoX;
	my $inter_rate_AtoX;

	my $blocks;
	my $blocks_AtoA;
	my $blocks_XtoX;
	my $blocks_AtoX;

	while(<SV>){
		chomp;
		$rates = $_;
#print "$line[0]\n";
#print "$rates\n";
		@rates = split/\s+/, $rates;
#print "*** @rates\n";
	}
	close SV;

	print "Rearrangement rates (off-diagnal / total syntenic blocks): $line[0] $rates[0]\n";
	print "Rearrangement rates A to A: $line[0] $rates[1]\n";
	print "Rearrangement rates X to X: $line[0] $rates[2]\n";
	print "Rearrangement rates A to X: $line[0] $rates[3]\n";

	print "Block length (Mb)(off-diagnal count / (off-diagnal length/1,000,000)): $line[0] $rates[4]\n";
	print "Block length A to A: $line[0] $rates[5]\n";
	print "Block length X to X: $line[0] $rates[6]\n";
	print "Block length A to X: $line[0] $rates[7]\n";

	print "Divergence time: $line[0] $time{$line[0]}\n";

	my $rate = $rates[0] / $time{$line[0]};
	my $rate_AtoA = $rates[1] / $time{$line[0]};
	my $rate_XtoX = $rates[2] / $time{$line[0]};
	my $rate_AtoX = $rates[3] / $time{$line[0]};
	my $block = $rates[4] / $time{$line[0]};
	my $block_AtoA = $rates[5] / $time{$line[0]};
	my $block_XtoX = $rates[6] / $time{$line[0]};
	my $block_AtoX = $rates[7] / $time{$line[0]};
	my $inter_rate = $rates[8] / $time{$line[0]};

	print OUT "$line[0]\t$time{$line[0]}\t$rates\t$rate\t$rate_AtoA\t$rate_XtoX\t$rate_AtoX\t$block\t$block_AtoA\t$block_XtoX\t$block_AtoX\t$inter_rate\n";
}
close LS;

print "Output file: 
 $ARGV[2]\n";
