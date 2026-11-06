#!/usr/bin/perl -w
use strict;

die "perl $0 [sv.list] [pairwise_time.tsv] [output file]
 Rearrangement rates are calculated as the number of rearrangements / total alignment (syntenic block) length of sp1 / divergence time\n" unless @ARGV == 3;

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

print OUT "spp\ttime\ttotal_align_len_in_Mb\trearrangements_counts\trearrangements_len\trates\tlen_rates\n";

while(<LS>){
	chomp;
	my $line = $_;
	my @line = split/\//, $line;
	my @spp = split/\_/, $line[0];
	print "Species pairs: $spp[0] - $spp[1]\n";

	open (SV, "$line") or die "$line $!\n";
	my $header = <SV>;
	print $header;

	my $align_len = 0;
	my $rearrange_len = 0;
	my $rearrange_count = 0;
	while(<SV>){
		chomp;
		my $sv_line = $_;
		my @sv_line = split /\s+/, $sv_line;

		# use sp1 alignment length as reference
		$align_len += $sv_line[6];
		
		if($sv_line[9] =~ /normal/){ # ignore normal syntenic blocks
			
		}else{
			$rearrange_count ++;
			$rearrange_len += $sv_line[6];
		}		
	}
	close SV;

	$align_len = $align_len / 1000000; # alignment length in Mb
	$rearrange_len = $rearrange_len / 1000000; # rearrangement length in Mb
	
	print "Rearrangement counts: $rearrange_count\n";
	print "Divergence time: $line[0] $time{$line[0]}\n";
	print "Total length of alignment: $line[0] $align_len\n";

	my $rate = $rearrange_count / $align_len / $time{$line[0]};
	my $len_rate = $rearrange_len / $align_len / $time{$line[0]};

	print OUT "$line[0]\t$time{$line[0]}\t$align_len\t$rearrange_count\t$rearrange_len\t$rate\t$len_rate\n";
}
close LS;

print "Output file: 
 $ARGV[2]\n";
