#/usr/bin/perl
use strict;

die "perl $0 [output from step1] [output]\n" unless @ARGV == 2;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[1]") or die "$ARGV[1] $!\n";

while(my $line = <IN>){
	chomp $line;
	my @line = split(/\t/,$line);

	my $len1;
	my $len2;
	my $dir;

	$len1 = abs($line[2] - $line[1]);
	$len2 = abs($line[5] - $line[4]);

	# orient the coordinates so spp1 is always from small to large
	if($line[1] < $line[2]){
		if($line[4] < $line[5]){
			$dir = "1";
			print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$len1\t$len2\t$dir\n";
		}else{
			$dir = "-1";
			print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[5]\t$line[4]\t$len1\t$len2\t$dir\n";
		}
	}elsif($line[1] > $line[2]){
		if($line[4] < $line[5]){
			$dir = "-1";
			print OUT "$line[0]\t$line[2]\t$line[1]\t$line[3]\t$line[4]\t$line[5]\t$len1\t$len2\t$dir\n";
		}else{
			$dir = "1";
			print OUT "$line[0]\t$line[2]\t$line[1]\t$line[3]\t$line[5]\t$line[4]\t$len1\t$len2\t$dir\n";
		}
	}
}
close IN;
close OUT;

print "Output: $ARGV[1]\n";
