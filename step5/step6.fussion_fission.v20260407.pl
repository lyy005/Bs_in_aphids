#!/usr/bin/perl -w
use strict;

die "perl $0 [step5.inter output] [X.list] [output file]
Updates:
  v20260307 - Calculate updated matrices such as length of inter-chr rearrangements and length per rearrangement 
  v20251104 - Add information of fussions and fissions in autosomes and Xs
  v20251031 - Infer fussions and fissions
         The idea is to find unique pairs of chr. If more than one chr is found, then it infers a fussion or fission event\n" unless @ARGV == 3;

my %x;

# find X chr
open (XL, "$ARGV[1]") or die "$ARGV[1] $!\n";
while(<XL>){
	chomp;
	my $line = $_;
	my @line = split /\s+/, $line;
	my $x = $line[-1];
	my @x = split/\,/, $x;

	foreach my $k (sort @x){
		$x{$k} = 1;
	}
}
close XL;

# find the diagnal chr (chr that has the highest synteny of each query chr to the ref chr)
# sp2 is the sp with more chromosomes
my %diagnal_len;
my %diagnal_count;
my %diagnal_chr;
my %offdiagnal_len;
my %offdiagnal_count;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";

my $header = <IN>;
print "$header\n";

# chr_pairs       chr_pair_syn_num        chr_pair_syn_len        ref_chr ref_num num_percent     ref_len len_percent
# ap1_ac1 9       3367990 ap1     106     0.0849056603773585      56058795        0.0600796003553055
# ap1_ac10        17      9067444 ap1     106     0.160377358490566       56058795        0.161748821036913

while(<IN>){
	chomp;
	my $line = $_;
	my @line = split /\s+/, $line;
	my $chr_pair = $line[0];
	my @chr_pair = split/\_/, $chr_pair;

	if($diagnal_len{$chr_pair[0]}){
		if($line[2] >= $diagnal_len{$chr_pair[0]}){
			$diagnal_len{$chr_pair[0]} = $line[2];
			$diagnal_count{$chr_pair[0]} = $line[1];
			$diagnal_chr{$chr_pair[0]} = $line[0];
			
		}
	}else{
		$diagnal_len{$chr_pair[0]} = $line[2]; # ap1 -> 9067444bp
		$diagnal_count{$chr_pair[0]} = $line[1]; # ap1 -> 17 rearrangements
		$diagnal_chr{$chr_pair[0]} = $line[0]; # ap1 -> ap1_ac10
	}
}
close IN;

print "Query-ref chr pairs with the longest syntenyic blocks:";
print "chr\tquery-ref\tnum_of_syntenic_blocks\tlength_of_blocks\n";
foreach my $k (sort keys %diagnal_chr){
	print "$k\t$diagnal_chr{$k}\t$diagnal_count{$k}\t$diagnal_len{$k}\n";
	
}

# by number
my $count = 0;
my $count_pairs = 0;
my $count_pairs_auto = 0;
my $count_pairs_x = 0;
my $count_pairs_autox = 0;

# by length
my $length = 0;
my $length_pairs = 0;
my $length_pairs_auto = 0;
my $length_pairs_x = 0;
my $length_pairs_autox = 0;

my (%sp1, %sp2);
my ($diagnal_len, $diagnal_count);
my ($diagnal_len_aa, $diagnal_count_aa);
my ($diagnal_len_xx, $diagnal_count_xx);
my ($diagnal_len_ax, $diagnal_count_ax);

my $offdiagnal_len = 0;
my $offdiagnal_count = 0;
my $offdiagnal_len_aa = 0;
my $offdiagnal_count_aa = 0;
my $offdiagnal_len_xx = 0;
my $offdiagnal_count_xx = 0;
my $offdiagnal_len_ax = 0;
my $offdiagnal_count_ax = 0;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OT1, ">$ARGV[2]") or die "$ARGV[2] $!\n";
open (OT2, ">$ARGV[2].log") or die "$ARGV[2].log $!\n";
$header = <IN>;
#chomp $header;
print "$header\n";

while(<IN>){
	chomp;
	$count_pairs ++;
	my $line = $_;
	my @line = split /\s+/, $line;

	my $chr_pair = $line[0];
	my @chr_pair = split/\_/, $chr_pair;

	# assign the chr to x or autosomes
	my ($chr1_loc, $chr2_loc);
	if(exists $x{$chr_pair[0]}){
		print "$chr_pair[0] is the X\n";
		$chr1_loc = "X";
	}else{
		$chr1_loc = "autosome";
	}

	if(exists $x{$chr_pair[1]}){
		print "$chr_pair[1] is the X\n";
		$chr2_loc = "X";
	}else{
		$chr2_loc = "autosome";
	}

	print OT2 "$chr_pair[0]\t$chr_pair[1]\t$chr1_loc\t$chr2_loc\n";

	# if the chr pair is the diagnal one / or the best match
	if($diagnal_chr{$chr_pair[0]} eq $chr_pair){
		$diagnal_len += $line[2];
		$diagnal_count += $line[1];
		
		# count auto-auto, x-x and auto-x rearrangements
		if(($chr1_loc eq "autosome")&&($chr2_loc eq "autosome")){
			$count_pairs_auto ++;

			$diagnal_len_aa += $line[2];
			$diagnal_count_aa += $line[1];

		}elsif(($chr1_loc eq "X")&&($chr2_loc eq "X")){
			$count_pairs_x ++;

			$diagnal_len_xx += $line[2];
			$diagnal_count_xx += $line[1];
		}else{
			$count_pairs_autox ++;

			$diagnal_len_ax += $line[2];
			$diagnal_count_ax += $line[1];
		}

	}else{
		$offdiagnal_len += $line[2];
		$offdiagnal_count += $line[1];

		# count auto-auto, x-x and auto-x rearrangements
		if(($chr1_loc eq "autosome")&&($chr2_loc eq "autosome")){
			$count_pairs_auto ++;

			$offdiagnal_len_aa += $line[2];
			$offdiagnal_count_aa += $line[1];

		}elsif(($chr1_loc eq "X")&&($chr2_loc eq "X")){
			$count_pairs_x ++;
	
			$offdiagnal_len_xx += $line[2];
			$offdiagnal_count_xx += $line[1];
		
		}else{
			$count_pairs_autox ++;

			$offdiagnal_len_ax += $line[2];
			$offdiagnal_count_ax += $line[1];
		}
	}

	# use sp1 as the ref
	if($sp1{$chr_pair[0]}){
		$sp1{$chr_pair[0]} .= ",".$chr_pair[1];
		$count ++;
	}else{
		$sp1{$chr_pair[0]} = $chr_pair[1];
	}

	# use sp1 as the ref
	if($sp2{$chr_pair[1]}){
		$sp2{$chr_pair[1]} .= ",".$chr_pair[0];
	}else{
		$sp2{$chr_pair[1]} = $chr_pair[0];
	}
}
close IN;

print OT2 "At least $count fussions and fissions\n";
print OT2 "Sp1 as ref\n";
my $sp1_chr = 0;
my $sp2_chr = 0;
foreach my $k (sort keys %sp1){
	print OT2 "$k\t$sp1{$k}\n";
	$sp1_chr ++;
}

print OT2 "Sp2 as ref\n";
foreach my $k (sort keys %sp2){
	print OT2 "$k\t$sp2{$k}\n";
	$sp2_chr ++;
}

print "Output: 
 At least $count fussions and fissions
 Details: $ARGV[2]\n";

my $inter_density = $count_pairs / ($sp1_chr * $sp2_chr);

my $inter_density_auto = $count_pairs_auto / ($sp1_chr * $sp2_chr);
my $inter_density_x = $count_pairs_x / ($sp1_chr * $sp2_chr);
my $inter_density_autox = $count_pairs_autox / ($sp1_chr * $sp2_chr);

print "Fusion and fission density: 
Number of uniq chr pairs / (#sp1 chr * #sp2 chr) = $count_pairs / ($sp1_chr * $sp2_chr) = $inter_density
Number of auto-auto pairs / (#sp1 chr * #sp2 chr) = $count_pairs_auto / ($sp1_chr * $sp2_chr) = $inter_density_auto
Number of X-X pairs / (#sp1 chr * #sp2 chr) = $count_pairs_x / ($sp1_chr * $sp2_chr) = $inter_density_x
Number of auto-X pairs / (#sp1 chr * #sp2 chr) = $count_pairs_autox / ($sp1_chr * $sp2_chr) = $inter_density_autox\n";

print OT2 "Fusion and fission density: 
Number of uniq chr pairs / (#sp1 chr * #sp2 chr) = $count_pairs / ($sp1_chr * $sp2_chr) = $inter_density
Density: $inter_density

Number of auto-auto pairs / (#sp1 chr * #sp2 chr) = $count_pairs_auto / ($sp1_chr * $sp2_chr) = $inter_density_auto
Auto-auto density: $inter_density_auto

Number of X-X pairs / (#sp1 chr * #sp2 chr) = $count_pairs_x / ($sp1_chr * $sp2_chr) = $inter_density_x
X-X density: $inter_density_x

Number of auto-X pairs / (#sp1 chr * #sp2 chr) = $count_pairs_autox / ($sp1_chr * $sp2_chr) = $inter_density_autox
Auto-X density: $inter_density_autox
";

#my $total_len;
#my $offdiagnal_len;
#my $diagnal_len;

#foreach my $k (sort keys %sp2_diagnal_chr){
#	$diagnal_len += $sp2_diagnal_len{$k};
#	$diagnal_count += $sp2_diagnal_count{$k};
#	print "$k\t$sp2_diagnal_chr{$k}\n";
#}

#foreach my $k (sort keys %sp2_offdiagnal_chr){
#	$offdiagnal_len += $sp2_offdiagnal_len{$k};
#	$offdiagnal_count += $sp2_offdiagnal_count{$k};
#	print "$k\t$sp2_offdiagnal_chr{$k}\n";
#}
my $offtarget_rate = ($offdiagnal_len / ($offdiagnal_len + $diagnal_len));
my $offtarget_rate_aa = ($offdiagnal_len_aa / ($offdiagnal_len + $diagnal_len));
my $offtarget_rate_xx = ($offdiagnal_len_xx / ($offdiagnal_len + $diagnal_len));
my $offtarget_rate_ax = ($offdiagnal_len_ax / ($offdiagnal_len + $diagnal_len));

my $offtarget_len;
if($offdiagnal_len == 0){
	$offtarget_len = 0;
}else{
	$offtarget_len = $offdiagnal_count / ($offdiagnal_len/1000000);
}

my $offtarget_len_aa;
if($offdiagnal_len_aa == 0){
	$offtarget_len_aa = 0;
}else{
	$offtarget_len_aa = $offdiagnal_count_aa / ($offdiagnal_len_aa/1000000);
}

my $offtarget_len_xx;
if($offdiagnal_len_xx == 0){
	$offtarget_len_xx = 0;
}else{
	$offtarget_len_xx = $offdiagnal_count_xx / ($offdiagnal_len_xx/1000000);
}

my $offtarget_len_ax;
if($offdiagnal_len_ax == 0){
	$offtarget_len_ax = 0;
}else{
	$offtarget_len_ax = $offdiagnal_count_ax / ($offdiagnal_len_ax/1000000);
}

print "
All chromosomes:
Inter-chromosomal rearrangement rate: 
Off-target length / Total syntenic length = $offdiagnal_len / ($offdiagnal_len + $diagnal_len) = $offtarget_rate
Off-target count / Off-target length = $offdiagnal_count / $offdiagnal_len = $offtarget_len Mb

A to A: 
Off-target length (A to A) / Total syntenic length = $offdiagnal_len_aa / ($offdiagnal_len + $diagnal_len) = $offtarget_rate_aa
Off-target count (A to A) / Off-target length (A to A) = $offdiagnal_count_aa / $offdiagnal_len_aa = $offtarget_len_aa Mb

X to X: 
Off-target length (X to X) / Total syntenic length = $offdiagnal_len_xx / ($offdiagnal_len + $diagnal_len) = $offtarget_rate_xx
Off-target count (X to X) / Off-target length (X to X) = $offdiagnal_count_xx / $offdiagnal_len_xx = $offtarget_len_xx Mb

A to X:
Off-target length (A to X) / Total syntenic length = $offdiagnal_len_ax / ($offdiagnal_len + $diagnal_len) = $offtarget_rate_ax
Off-target count (A to X) / Off-target length (A to X) = $offdiagnal_count_ax / $offdiagnal_len_ax = $offtarget_len_ax Mb
";


print OT1 "Interchr_block_rate\tInterchr_block_rate_AtoA\tInterchr_block_rate_XtoX\tInterchr_block_rate_AtoX\tAvg_number_of_breaks_perMb\tAvg_number_of_breaks_perMb_AtoA\tAvg_number_of_breaks_perMb_XtoX\tAvg_number_of_breaks_perMb_AtoX
$offtarget_rate\t$offtarget_rate_aa\t$offtarget_rate_xx\t$offtarget_rate_ax\t$offtarget_len\t$offtarget_len_aa\t$offtarget_len_xx\t$offtarget_len_ax\n";
