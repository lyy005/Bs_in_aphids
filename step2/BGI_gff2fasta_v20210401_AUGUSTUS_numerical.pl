#!/usr/bin/perl

use strict;

die "perl $0 [gff file] [fasta file] [cds output]\nAuthor: Yiyuan LI, liyiyuan\@genomics.org.cn
 For Hormaphis cornu braker annotation
 Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" unless (@ARGV == 3);

#my $pos_file = shift;
#my $seq_file = shift;

my %gene;

read_gff($ARGV[0],\%gene);

open IN, "$ARGV[1]" or die "$!\n";
open CDS, ">$ARGV[2]" or die "$!\n";
open CDSLOG, ">$ARGV[2].log" or die "$!\n";

$/=">"; <IN>; $/="\n";	
while (<IN>) {
	my $output;
	my $chr=$1 if(/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$/="\n";

	my $seq_len=length($seq);
	my $chr_pp=$gene{$chr};
	#foreach my $gene (sort {($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0]} keys %$chr_pp) {
	# sort by gene IDs (g1 g2 etc)
	foreach my $gene (sort {($a =~ /g(\d+)/)[0] <=> ($b =~ /g(\d+)/)[0]} keys %$chr_pp) {
		my @strand=@{$chr_pp->{$gene}{strand}};
		next if(!exists $chr_pp->{$gene}{exon});
		my @exon = @{$chr_pp->{$gene}{exon}};
		
		my $mrna;
#		my ($left_leng, $right_leng);
#		$left_leng = $exon[0][0]-1 if($left_leng > $exon[0][0]-1);
#		$right_leng = $seq_len - $exon[-1][1] if($right_leng > $seq_len - $exon[-1][1]);
#		$mrna .= substr($seq,$exon[0][0]-$left_leng-1,$left_leng) if($left_leng);
		for (my $i=0; $i<@exon; $i++) {
#			$mrna .= substr($seq,$exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1);	
			my $cds = substr($seq,$exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1);
			print CDSLOG "$gene $i $cds\n"; 
			#$cds = Complement_Reverse($cds) if($strand[$i] eq '-'); # modification 1
			$mrna .= $cds;
		}
		$mrna = Complement_Reverse($mrna) if($strand["0"] eq '-');  # modification 1
		print CDSLOG "$gene combined $mrna\n"; 
#		$mrna .= substr($seq,$exon[-1][1],$right_leng); # if($right_leng);
#		$mrna = Complement_Reverse($mrna) if($strand eq '-');
		Display_seq(\$mrna);
		my $mark = "$gene  [CDS]  locus=$chr:$exon[0][0]:$exon[-1][1]";
		$output .= ">".$mark."\n".$mrna;
		
	}
	print CDS $output;
}

close(IN);




####################################################
################### Sub Routines ###################
####################################################

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################


#############################################
sub Complement_Reverse{
	my $seq=shift;
	$seq=~tr/AGCTagct/TCGAtcga/;
	$seq=reverse($seq);
	return $seq;

}
#############################################


sub read_gff{
	my $file=shift;
	my $ref=shift;
	open (IN,$file) || die ("fail open $file\n");
	open LOG, ">$file.log" or die "$!\n";

	while (<IN>) {
		next if(/^\#/);
		s/^\s+//;
		s/\s+$//;
		my @t = split(/\t/);
		my $tname = $t[0];
		my ($qname, $gene_name);
#		if ($t[2] eq 'mRNA' || $t[2] eq 'CDS') {
		if ($t[2] eq 'CDS') {
			$qname = $1 if($t[8] =~ /Parent=([^\;]+)/);
			$gene_name = $1 if($t[8] =~ /Dbxref=GeneID\:([^\,]+)\,/);
			#$qname = $gene_name."-".$qname;
			$qname = $qname;

			my $id = $1 if($t[8] =~ /ID=([^\;]+)\;/);

#			$ref->{$tname}{$qname}{strand} = $t[6];
			push @{$ref->{$tname}{$qname}{strand}}, $t[6];
			push @{$ref->{$tname}{$qname}{exon}}, [$t[3],$t[4]];

			print LOG "$tname\t$id\t$qname\t$t[3]\t$t[4]\t$t[6]\n";
		}
	}
	close IN;
	close LOG;
	
}
