# Bs_in_aphids
Pipeline for identifying supernumerary B chromosomes in sumac gall aphids

## Outline: 
- [1 - Improving the assembly of Bs](#1---improving-the-assembly-of-bs)
 - [1.1 - Improving the assembly of Bs](#1---improving-the-assembly-of-bs)

## 1 - Genome assembly and annotation
To simplify the B chromosome assembly process, we first identified autosomes and Xs based on sequencing homology. Then we assembled the B-derived reads to improve the assembly of Bs.

### 1.1 - Identifying autosomes and Xs based on the genome of *Schlechtendalia chinensis*
We first aligned hifiasm results to the top 13 scaffolds (chromosomes) of *Schlechtendalia chinensis* using LAST. Then we picked the top scaffolds that matched the chromosomes of *Schlechtendalia chinensis*
```
conda activate last_v1648

# ref: https://github.com/aaannaw/Myospalax-CFi/blob/main/CFi_identification/Ancestral_karyotype_reconstruction.sh
# ref: https://gitlab.com/mcfrith/last/-/blob/1648/doc/last-cookbook.rst

# human and chimp example (for example, s5h1 vs s1 and s5h1 vs s10h1 )
#lastdb -P10 -c -U400 -uRY128 ref_db ref.fa
#last-train -P10 --revsym -C2 ref_db q.fa > hc.train
#lastal -P10 -D1e9 -C2 --split-f=MAF+ -p hc.train ref_db q.fa | last-split -r | maf-linked - > aln.maf

# Aligning more distantly related genomes (for example, s5h1 vs s6h1)
lastdb -P10 -c ref_db ref.fa
last-train -P10 --revsym -C2 ref_db q.fa > hl.train
lastal -P10 -m2 -D1e9 -C2 --split-f=MAF+ -p hl.train ref_db q.fa | last-split -r | maf-linked - > aln.maf

conda activate last_v1542
last-dotplot aln.maf aln.maf.png --rot1=v --rot2=h
```

The resulting dotplot looks like this:
<img width="527" height="1000" alt="aln maf" src="https://github.com/user-attachments/assets/d03513a8-0e4c-46ca-a581-951b75c86fa1" />


#### Identify Xs based on resequencing data
To further confirm the identified Xs, we sequenced males and females of *Kaburagia rhusicola rhusicola*, 

#### Map HiFi reads to hifiasm scaffolds
Using samtools version 1.21 and minimap2 version 2.28
```
conda activate samtools_v1.21

# cat h1.rescaf.fa h2.rescaf.fa > genome.fa
minimap2 -a genome.fa -o genome.sam -t 20 -x map-hifi hifi.filt.fastq.gz
samtools view -h -b -S genome.sam -o genome.bam --threads 20
samtools sort genome.bam -o genome.sorted.bam --threads 20
samtools index genome.sorted.bam
```

Use SAMtools to keep reads that did not align to autosomes and Xs
```
conda activate samtools_v1.21

samtools view -bS genome.sorted.bam -o B_reads.bam h1tg000021l h1tg000001l h1tg000010l h2tg000015l h1tg000013l h2tg000025l h1tg000026l h2tg000008l h1tg000015l h1tg000017l h2tg000017l h1tg000022l h1tg000030l h1tg000003l h1tg000023l h1tg000027l h1tg000024l h1tg000028l h1tg000042l h1tg000048l h1tg000029l h1tg000054l h1tg000039l h1tg000034l h1tg000055l h1tg000059l h1tg000056l h1tg000060l h1tg000076l h1tg000071l h1tg000065l h1tg000050l h1tg000049l h1tg000057l h1tg000025l h1tg000041l h1tg000052l h1tg000036l h1tg000078l h1tg000043l h1tg000031l h1tg000074l h1tg000067l h1tg000037l h1tg000072l h1tg000062l h1tg000044l h1tg000068l h1tg000069l h1tg000051l h1tg000064l h1tg000035l h1tg000075l h1tg000033l h1tg000073l h1tg000077l h1tg000061l h1tg000045l h1tg000040l h1tg000038l h1tg000063l h1tg000070l h1tg000053l h1tg000046l h1tg000066l h1tg000047l h1tg000058l h1tg000032l h2tg000020l h2tg000024l h2tg000023l h2tg000027l h2tg000030l h2tg000022l h2tg000032l h2tg000018l h2tg000031l h2tg000036l h2tg000029l h2tg000034l h2tg000035l h2tg000039l h2tg000042l h2tg000037l h2tg000033l h2tg000028l h2tg000040l h2tg000045l h2tg000046l h2tg000021l h2tg000044l h2tg000026l h2tg000043l h2tg000041l h2tg000038l --threads 20

samtools fastq B_reads.bam | gzip B_reads.fastq.gz
```

#### Assemble B-derived reads using Verkko
```
conda activate verkko_v2.2.1
verkko -d run1_asm --hifi hifi.filt.fastq.gz --telomere-motif CCTAA
```

## 2 - Assessing the completeness of the genome and Bs
To assess the completeness of the genome and the Bs, we checked the completeness of telomeres. We used two approaches: 

### 2.1 - seqtk
seqtk searches telomeres are the very end of each scaffold:
```
# search for telomeres at ends of each scaffolds
conda activate seqtk_v1.4-r122
seqtk telo -m CCTAA genome.fa > genome.fa.telo

# calculate the length of the assembly
perl ./step1_assembly/length_distribution_V1.5.pl genome.fa

# summarize the location of telomeres to a table
perl ./step1_assembly/summarize_telomere.v0.92.pl genome.fa.telo genome.fa.length genome.fa.telo.summary
```
### 2.2 - TRF
On top of telomeres found by seqtk, TRF also searches for telomeres in the middle or close to the end of each scaffold:
```
conda activate trf_v4.10.0
trf genome.fa 2 7 7 80 10 50 500 -f -d -m
```

## 3 - Characterize Bs in aphids
We also characterized other features of the genomes and the Bs

### 3.1 - Annotation of ribosomal RNAs
We used barrnap (https://github.com/tseemann/barrnap) to annotated rRNAs. 
```
conda activate barrnap
barrnap --quiet s5h1.fa --kingdom euk
```

### 3.2 - Annotation of repeats
To annotate repetitive content, we used EarlGrey. 

We first annotated each genome separately:
```
conda activate earlgrey_v6.3.2
earlGrey -g kaburagiaRhusicola.fa -s kaburagiaRhusicola -r arthropoda -o output -t 20 -c yes -e yes -d yes
```

Then we combined all the .strained outputs:
```
for i in s10h1 s10h2 s11h1 s11h2 s1 s3h1 s3h2 s4 s5h1 s5h2 s6
do 
 echo $i
 perl add_spp_name.pl $i\.strained $i $i\.renamed.strained
done
```
Clustered the strained sequences (ref: https://github.com/TobyBaril/EarlGrey/issues/216):

```
cat *.renamed.strained > combined.renamed.strained

conda activate mmseqs2_v18
mmseqs easy-cluster combined.renamed.strained combined.renamed.strained.easyCluster tmp --min-seq-id 0.8 -c 0.8 --cov-mode 1 -v 2 --threads 10 --cluster-reassign
```

Lastly, re-annotated all the genomes using the same library:
```
conda activate earlgrey_v6.3.2

earlGreyAnnotationOnly -g kaburagiaRhusicola.fa -s kaburagiaRhusicola -o output -t 20 -l combined.renamed.strained.easyCluster_rep_seq.fasta
```

### 3.3 - DNA Methylation Analysis
To identify DNA methylcations, we used pbmm2 (https://github.com/PacificBiosciences/pbmm2) and aligned_bam_to_cpg_scores version 3.0 (https://github.com/PacificBiosciences/pb-CpG-tools). 

```
conda activate pbmm2_1.17.0

pbmm2 index genome.fa genome.mmi
pbmm2 align genome.mmi hifi_reads.bam aligned.bam --preset CCS --sort --log-level DEBUG --num-threads 20

aligned_bam_to_cpg_scores --bam aligned.bam --output-prefix cpg_denovo --threads 20 --modsites-mode denovo --ref genome.fa
```

### 3.4 - TAD

```
step0 - mapping
conda activate hicexplorer_v3.7.6

bwa mem -A 1 -B 4 -E 50 -L 0 -t 40 genome.fa mate_R1.fastq.gz 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam
bwa mem -A 1 -B 4 -E 50 -L 0 -t 40 genome.fa mate_R2.fastq.gz 2>>mate_R2.log | samtools view -Shb - > mate_R2.bam
hicBuildMatrix --samFiles mate_R1.bam mate_R2.bam \
               --binSize 10000 \
               --restrictionSequence GATC \
               --danglingSequence GATC \
               --restrictionCutFile cut_sites.bed \
               --threads 40 \
               --inputBufferSize 100000 \
               --outBam ./hicMatrix/hic.bam \
               -o ./hicMatrix/hic_matrix.h5 \
               --QCfolder ./hicMatrix/hicQC \
               --minMappingQuality 15

# step1 - correction
hicCorrectMatrix diagnostic_plot -m hic_matrix.h5 --perchr -o hic_matrix.diagnostic_plot.pdf

hicCorrectMatrix correct -m hic_matrix.h5 -o hic_KRcorrected.h5 --correctionMethod KR
hicCorrectMatrix diagnostic_plot -m hic_KRcorrected.h5 --perchr -o hic_KRcorrected.diagnostic_plot.pdf

## change resolution from 10kb to 30kb
hicMergeMatrixBins --matrix hic_KRcorrected.h5 --outFileName hic_KRcorrected_30kb.h5 --numBins 3

# step2 - cis vs trans
hicConvertFormat --matrices hic_KRcorrected_30kb.h5 --inputFormat h5 --outputFormat ginteractions --outFileName hic_KRcorrected_30kb.ginteractions

## filter interactions with counts >= 3
less -S hic_KRcorrected_30kb.ginteractions.tsv | awk '$7 >= 3' | cut -f 1,4 | sort | uniq -c > hic_KRcorrected_30kb.ginteractions.counts

for i in B1 B2 B3 B4 B5
do
 perl step3_add_length.pl hic_KRcorrected_30kb.ginteractions.counts genome.fa.length $i hic_KRcorrected_30kb.ginteractions.$i\.counts
done

```



