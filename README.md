# Bs_in_aphids
Pipeline for identifying supernumerary B chromosomes in sumac gall aphids

Photos of the sumac galls induced by these aphids. Scale bar represents 1 cm. 
<img width="2718" height="850" alt="image" src="https://github.com/user-attachments/assets/d93c75b2-099a-46a9-972b-76072f1c7e57" />



## Outline: 
- [1 - Mitochondrial genome assembly and species identification](#1---mitochondrial-genome-assembly-and-species-identification)
  - [1.1 - Quality control](#11---quality-control)
  - [1.2 - Assembly of mitochondrial genomes from Illumina data using MitoZ](#12---assembly-of-mitochondrial-genomes-from-illumina-data-using-mitoz)
  - [1.3 - Annotation of mitochondrial genomes from NCBI using MitoZ](#13---annotation-of-mitochondrial-genomes-from-ncbi-using-mitoz)
  - [1.4 - Phylogeny reconstruction using mitochondrial genes](#14---phylogeny-reconstruction-using-mitochondrial-genes)
- [2 - Genome assembly and annotation](#2---genome-assembly-and-annotation)
  - [2.1 - Generating CCS reads](#21---generating-ccs-reads)
  - [2.2 - Genome assembly using hifiasm](#22---genome-assembly-using-hifiasm)
  - [2.3 - Scaffolding using HiC data](#23---scaffolding-using-hic-data)
  - [2.4 - Manual curation](#24---manual-curation)
  - [2.5 - BUSCO assessment](#25---busco-assessment)
  - [2.6 - Gene annotation](#26---gene-annotation)
  - [2.7 - Telomere annotation](#27---telomere-annotation)
  - [2.8 - Ribosomal DNA annotation](#28-ribosomal-dna-annotation)
  - [2.9 - Transposable element annotation](#29---transposable-element-annotation)
- [3 - X chromosome assignment](#3---x-chromosome-assignment)
  - [3.1 - Identifying Xs of *K. r. ovogallis* based on the sequencing depth of males and females](#31---assigning-x-chromosomes-of-k-r-ovogallis-based-on-the-sequencing-depth-of-males-and-females)
  - [3.2 - Identifying Xs of other species based on the genome synteny with *K. r. ovogallis*](#32---assigning-x-chromosomes-of-other-species-based-on-the-genome-synteny-with-k-r-ovogallis)
- [4 - Phylogeny reconstruction](#4---phylogeny-reconstruction)
  - [4.1 - OrthoFinder](#41---ortholog-inference)
  - [4.2 - Phylogenetic tree reconstruction](#42---phylogenetic-tree-construction)
  - [4.3 - Divergence time estimation](#43---divergence-time-estimation)
- [5 - Chromosomal rearrangement rate estimation](#5---chromosomal-rearrangement-rate-estimation)
  - [5.1 - Divergence time between species pairs](#51---divergence-time-between-species-pairs)
  - [5.2 - Estimating intra-chromosomal rearrangement rates](#52---estimating-intra-chromosomal-rearrangement-rates)
  - [5.3 - Estimating inter-chromosomal rearrangement rates](#53---estimating-inter-chromosomal-rearrangement-rates)
  - [5.4 - Visualization of rates](#54---visualization-of-rates)
- [6 - B chromosome assembly]
- [7 - Estimation of *dN/dS* and gene expression](#7---dnds-and-gene-expression-of-b-linked-genes)
  - [7.1 - *dN/dS* ](#71---measuring-dnds)
  - [7.2 - TPM](#72---estimating-tpm)
- [8 - 3D organization of B chromosomes](#8---3d-organization-of-b-chromosomes)
  - [8.1 - Topologically associating domain](#81---tad)
  - [8.2 - Trans-interactions between chromosomes](#82---trans-interactions)
- [9 - Methylation of B chromosomes](#9---methylation-in-b-chromosomes)
  - [9.1 - Methylation detection](#91---methylation-detection)
  - [9.2 - CpGs detected in the genome](92---cpg-sites-in-the-genome)
  - [9.3 - Exon, intron, and intergenic regions](#93---identifying-exon-intron-and-intergenic-regions)
  - [9.4 - Summarization and visualization of CpGs]()

## 1 - Mitochondrial genome assembly and species identification
### 1.1 - Quality control
We first removed low-quality bases and sequencing adapters using Trimmomatic. The combined_MGI.fa file combines Illumina and MGI sequencing adapters. We used the same parameters for quality control of other sequencing libraries. 
```
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 combined.1.fq.gz combined.2.fq.gz \
          pe.1.fq.gz se.1.fq.gz pe.2.fq.gz se.2.fq.gz \
          ILLUMINACLIP:./step_1/combined_MGI.fa:2:30:10:8:TRUE \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 10
```

### 1.2 - Assembly of mitochondrial genomes from Illumina data using MitoZ
We assembled clean pair-end reads using MitoZ version 3.6 (https://github.com/linzhi2013/MitoZ). For the aphid mitochondrial genomes, we used the invertebrate mitochondrial genetic code (“--clade Arthropoda --requiring_taxa Arthropoda” and metaSPAdes (“--assembler spades”). 

```
singularity run ~/bin/mitoZ/MitoZ_v3.6.sif mitoz all --fq1 pe.1.fq.gz --fq2 pe.2.fq.gz \
            --outprefix mitoz --thread_number 96 --workdir ./ --clade Arthropoda \
            --requiring_taxa Arthropoda --assembler spades --memory 1400
```

### 1.3 - Annotation of mitochondrial genomes from NCBI using MitoZ
For other mitochondrial genomes from NCBI, we ran MitoZ to annotate the protein-coding genes so that the formats of mitochondrial genes were consistent for the mitochondrial genomes from our study and from NCBI. 

```
singularity run ~/bin/mitoZ/MitoZ_v3.6.sif mitoz annotate --workdir ./ --outprefix NCBI \
          --thread_number 20 --fastafiles NCBI.fasta --species_name NCBI --clade Arthropoda
```

### 1.4 - Phylogeny reconstruction using mitochondrial genes
After mitoZ assembly and annotation, all mitochondrial genes should look as follows. Each gene name is species name + ";" + gene name + ";" + length + coordinates on the genome
```
>Bapi;COX1;len=1536;[0:1536](+)
ATAAAATGAATTTTTTCAACTAATCATAAAAATATTGGAACTTTATATTTTTTATTTGGTTTATGATCTGGAATAATTGGATCATCACTTAGAATATTA
>Bapi;COX2;len=672;[1601:2273](+)
ATATCTTGATTAAAACTAAGATTTCAAAATAGAAATTCACCATTAATAGAACAATTAATCTTCTTTCATGATCATACTATCTTTATTATTATCATAATT
>Bapi;ATP8;len=159;[2419:2578](+)
ATACCACAAATAGCTCCTATAAATTGATTAATAATATT
```

We then combined this file from all the species into combined.fasta. We first removed ":" characters. Then use a Perl script to separate sequences by genes:
```
less -S ./step1/combined.fasta | perl -ne 's/\:\>/\:/g; print; ' > cds.noCaret.fasta

perl ./step1/step1.split_by_gene.pl cds.noCaret.fasta ./cds
```

Then we aligned each gene based on codon sequences (https://github.com/lyy005/codon_alignment) and removed poorly aligned sites using Gblocks version 0.91b. 
```
for i in ATP6 ATP8 COX1 COX2 COX3 CYTB ND1 ND2 ND3 ND4 ND4L ND5 ND6
do
 echo $i
 # codon alignment
 perl ~/bin/codon_alignment_v1.5/codon_alignment.pl $i\.nt.fasta 5 0 > $i\.aln.log

 # Gblocks
 Gblocks $i\.nt.fasta.nt.cleanup.aln -t=c -b4=10 -b5=n -e=-gb > $i\.GBlock.log

 # convert gblocks output to fasta formats
 perl ./step1/step2.gblock_convert.pl $i\.nt.fasta.nt.cleanup.aln-gb > $i\.GBlocks.fas
done
```

Lastly, we built phylogenetic trees using IQ-Tree version 3.0.1. To do that, we need two files: 1. a file with all the species names, and 2. a file with the path to each gene alignment. 

The spp.list file is under ./step1/

```
# make the gene list file: 
ls ./cds/*.GBlocks.fas > gene.list

# concatenate genes together
perl step3.Concatenation_v1.0.pl gene.list concatenated.fasta spp.list

# make the tree
~/bin/iqtree/iqtree-3.0.1-Linux/bin/iqtree3 -s concatenated.fasta -B 1000 -T 10 -o Bapi --prefix iqtree3_output

```

## 2 - Genome assembly and annotation
To simplify the B chromosome assembly process, we first identified autosomes and Xs based on sequencing homology. Then we assembled the B-derived reads to improve the assembly of Bs.

### 2.1 - Generating CCS reads
*K. r. ovogallis*, *K. r. ovatirhusicola*, *K. r. rhusicola*, and *K. r. sp.* were sequenced on PacBio II platforms and raw subreads were generated. To retrieve the Circular Consensus Sequencing (CCS) reads, we used ccs (https://github.com/PacificBiosciences/ccs) to convert subreads to CCS reads. Reads from each pass were assembled with hifiasm. The pass with the highest N50 and the most complete BUSCOs was used for the final assembly. 

Other species were sequenced on PacBio II platforms. Therefore, CCS reads were directly generated and assembled using hifiasm. 


```
for i in 3 4 5 6 7 8 9 10
do
 ccs example.subreads.bam ccs$i\.fastq.gz --minPasses $i --reportFile ccs$i\_report.txt --num-threads 40
done
```

### 2.2 - Genome assembly using hifiasm
We used hifiasm version 0.25.0 to assemble the HiFi reads with dual scaffolding, aphid telomere sequences, and HiC data. 

Before assembly, HiC data were also filtered using Trimmomatic using the same command as the Quality control process in [1.1 - Quality control](#11---quality-control). 

```
~/bin/hifiasm_v0.25.0/hifiasm -o ccs --dual-scaf --telo-m CCTAA --primary -t 96 --h1 HiC_R1.fq.gz --h2 HiC_R2.fq.gz hifi.filt.fastq.gz

# converting .gfa files to .fasta files
awk '/^S/{print ">"$2;print $3}' ccs.hic.hap1.p_ctg.gfa > ccs.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ccs.hic.hap2.p_ctg.gfa > ccs.hic.hap2.p_ctg.fa
```

### 2.3 - Scaffolding using HiC data
The outputs of hifiasm were scaffolded based on HiC contacts using the Arima Genomics mapping pipeline (https://github.com/ArimaGenomics/mapping_pipeline) and YaHS version 1.2.2 (https://github.com/c-zhou/yahs). Each haplotype was scaffolded independently. draft.fa file was the output of hifiasm assembly (either ccs.hic.hap1.p_ctg.fa or ccs.hic.hap2.p_ctg.fa).

<details>
<summary> Click to show the HiC mapping code </summary>
	
```
# Step 0: Index reference
~/bin/bwa/bwa_v0.7.17/bwa-0.7.17/bwa index -a bwtsw -p draft.fa draft.fa

# Step 1.A: FASTQ to BAM (forward HiC reads)
~/bin/bwa/bwa_v0.7.17/bwa-0.7.17/bwa mem -t 4 draft.fa hic_R1.fastq.gz | samtools view -@ 4 -Sb - > hic_R1.bam

# Step 1.B: FASTQ to BAM (reverse HiC reads)
~/bin/bwa/bwa_v0.7.17/bwa-0.7.17/bwa mem -t 4 draft.fa hic_R2.fastq.gz | samtools view -@ 40 -Sb - > hic_R2.bam

# Step 2.A: Filter 5' end (forward HiC reads)
samtools view -h hic_R1.bam | perl ~/Arima_mapping_pipeline/filter_five_end.pl | samtools view -Sb - > hic_R1.filter.bam

# Step 2.B: Filter 5' end (2nd)
samtools view -h hic_R2.bam | perl ~/Arima_mapping_pipeline/filter_five_end.pl | samtools view -Sb - > hic_R2.filter.bam

# Step 3A: Pair reads & mapping quality filter
perl ~/Arima_mapping_pipeline/two_read_bam_combiner.pl hic_R1.filter.bam hic_R2.filter.bam samtools 10 | samtools view -bS -t draft.fa.fai - |samtools sort -@ 40 -o hic.combined.bam -

# Step 3B: Add read group
java -Xmx4G -Djava.io.tmpdir=temp/ -jar ~/bin/picard/picard_v2.27.3/picard.jar AddOrReplaceReadGroups INPUT=hic.combined.bam OUTPUT=hic.combined.readgroup.bam ID="hic" LB="hic" SM=aphids PL=ILLUMINA PU=none

# Step 4: Mark duplicates
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar ~/bin/picard/picard_v2.27.3/picard.jar MarkDuplicates INPUT=hic.combined.readgroup.bam OUTPUT=hic.combined.readgroup.dereplicate.bam METRICS_FILE=metrics.hic.combined.readgroup.dereplicate.txt TMP_DIR=temp/ ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index hic.combined.readgroup.dereplicate.bam

perl ~/Arima_mapping_pipeline/get_stats.pl hic.combined.readgroup.dereplicate.bam > hic.combined.readgroup.dereplicate.bam.stats

# scaffolding
samtools faidx draft.fa ~/bin/yahs/yahs_v1.2.2/yahs --no-contig-ec -e GATC draft.fa hic.combined.readgroup.dereplicate.bam
```

</details>

We used the final YaHS output, yahs.out_scaffolds_final.agp, as the reference for scaffolding. But we also performed manual curation to confirm that the corrections suggested by YaHS were correct. 

### 2.4 - Manual curation
We manually curated scaffolds using CurationPretext version 1.6.1 (https://github.com/sanger-tol/curationpretext) and PretextView version 1.0.3 (https://github.com/sanger-tol/PretextView).

We first converted HiC .fastq.gz files to .cram files using SAMtools and HiFi .fastq.gz files to .fasta format. We ran the CurationPretext pipeline on the scaffolds from both haplotypes 1 and 2 from hifiasm.

```
samtools import -@10 -r ID:aphid -r CN:aphid -r PU:hic -r SM:sample6 hic_R1.fastq.gz hic_R2.fastq.gz -o hic.cram
samtools index hic.cram

conda activate seqkit_v2.2.0
seqkit fq2fa hifi.filt.fastq.gz -o hifi.fasta.gz --threads 10

cat ccs.hic.hap1.p_ctg.fa ccs.hic.hap2.p_ctg.fa > genome.fa
```

Then we made a parameters file:
```
input: "~/PretextView_hifiasm/genome.fa"
outdir: "./output/"
teloseq: "TTAGG"
sample: "s5"
reads: ["~/PretextView_hifiasm/pacbio/hifi.fasta.gz"]
read_type: "hifi"
cram: ["~/PretextView_hifiasm/hic/hic.cram"]
```

We ran the CurationPretext pipeline:
```
nextflow run ~/bin/curationpretext/v1.6.1 -profile singularity -params-file params.yaml
```

### 2.5 - BUSCO assessment
We measured the quality of the genome assembly using BUSCOs (http://busco.ezlab.org/):
```
conda activate busco_v5.4.2

busco -i genome.fa -o genome_hemiptera -m genome -l hemiptera_odb10 -c 10 -f
```

### 2.6 - Gene annotation
To annotate genes, we followed the annotation pipeline of BRAKER version 3.0.8. 
We first soft-masked the genome using tantan version 51. Then we aligned RNA-seq data to the genome using HISAT2 version 2.2.1. The RNA-seq data were filtered the same way as [1.1](#11---quality-control). Then we used BRAKER to annotate the genome based on both RNA-seq evidence and protein evidence of related aphid species (proteins_NCBI_GCF_aphids_v20210211_noStopCodon.fasta file in the step2 folder).
```
tantan genome.fa -f 0 > genome.masked.fa

hisat2-build genome.masked.fa genome.masked.fa
hisat2 -x genome.masked.fa -1 pe.1.fq.gz -2 pe.2.fq.gz -U se.1.fq.gz,se.2.fq.gz -S combined.sam --phred33 -t --threads 20

singularity exec ~/bin/braker/BRAKER-3/braker3.sif braker.pl --genome=genome.masked.fa \
	--bam=combined.sorted.bam \
	--prot_seq=proteins_NCBI_GCF_aphids_v20210211_noStopCodon.fasta \
	--softmasking \
	--threads 20 \
	--workingdir=braker3_output \
	--species example \
	--gff3 \
	--busco_lineage=hemiptera_odb10
```

We then renamed the genes so their IDs increased in order along the chromosomes. 
<details>
<summary> Click to show the code to reformat the outputs from BRAKER </summary>
	
```
perl ./step2/length_distribution_V1.5.pl genome.masked.fa

# sort by scaffold
perl ./step2/sort_gtf_by_scaffold_ids.v0.9.pl braker.gtf genome.masked.fa.length braker.sorted.gtf
perl ./step2/rename_gtf_by_scaffold_order.pl braker.sorted.gtf braker.renamed.gtf


# GTF to GFF
conda activate braker_v20210330_mamba
gtf2gff.pl < braker.renamed.gtf --gff3 --out braker.renamed.gff3

# add scaf names to genes
less -S braker.renamed.gff3 | perl -ne '@a=split; s/ID\=/ID\=$a[0]\_/; s/Parent\=/Parent\=$a[0]\_/; print; ' > braker_hap_named.gff3

# gff to cds
perl ./step2/BGI_gff2fasta_v20210401_AUGUSTUS_numerical.pl braker_hap_named.gff3 genome.masked.fa braker.cds.fa

# cds to pep
perl ./step2/cds2aa_stopCodonInAsterisks.pl braker.cds.fa NUC braker.pep.fa

# remove stop codons
perl ./step2/remove_genes_with_stop_codons.pl braker.pep.fa braker.pep.noStop.fa

# longest isoform
perl ./step2/find_longest_protein_Braker_v20210403b.braker.pl braker.pep.noStop.fa braker.pep.noStop.fa.longest
```

</details>

### 2.7 - Telomere annotation
To assess the completeness of assembled chromosomes, we annotated telomeric sequences using two approaches: seqtk version 1.4 (https://github.com/lh3/seqtk) and telomere version 0.0.4 (https://github.com/sanger-tol/telomere). 

seqtk searches for telomeres at the very end of each scaffold:
```
# search for telomeres at ends of each scaffold
conda activate seqtk_v1.4-r122
seqtk telo -m CCTAA genome.fa > genome.fa.telo

# calculate the length of the assembly
perl ./step2/length_distribution_V1.5.pl genome.fa

# summarize the location of telomeres to a table
perl ./step2/summarize_telomere.v0.92.pl genome.fa.telo genome.fa.length genome.fa.telo.summary
```

telomere searches for the percentage of telomeres in sliding windows across the genome.
```
~/bin/telomere-0.0.4/find_telomere genome.fa TTAGG > genome.telomere
cat genome.telomere | awk '{print $1"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' - > genome.telomere.reform
sdust genome.fa > genome.sdust

java -cp ~/bin/telomere-0.0.4/telomere.jar SizeFasta genome.fa > genome.lens 

# FindTelomereWindows writes prefix.telomere.bed and prefix.windows
# (and with --split also fwd/rev beds + fwd/rev windows).
# Lowering threshold to 0.10 (10%) from the initial 0.40 (40%)
java -cp ~/bin/telomere-0.0.4/telomere.jar FindTelomereWindows genome.telomere 99.9 0.1
java -cp ~/bin/telomere-0.0.4/telomere.jar FindTelomereBreaks genome.lens genome.sdust genome.telomere > genome.breaks

# check the percentage of telomeres in each window
less -S genome.windows | sort -V -k1,1 -k2,2n | less -S
```

### 2.8 Ribosomal DNA annotation
Aphid ribosomal DNA arrays are located on the X chromosome (Blackman and Spence 1996). To locate the rDNA arrays in aphids with multiple Xs, we annotated the rDNAs using barrnap (https://github.com/tseemann/barrnap). 

```
conda activate barrnap
barrnap --quiet genome.fa --kingdom euk
```

### 2.9 - Transposable element annotation
To understand the distribution of transposable elements, we annotated transposable elements using Earl Grey version 6.3.0 (https://github.com/TobyBaril/EarlGrey). To ensure consistent annotation across genomes, we then combined and clustered the de novo curated elements using mmseqs version 18. This method was suggested by Earl Grey: https://github.com/TobyBaril/EarlGrey/issues/216

```
# First, we annotated each genome using Earl Grey
earlGrey -g kaburagiaRhusicola.fa -s kaburagiaRhusicola -r arthropoda -o output -t 20 -c yes -e yes -d yes

# Secondly, we renamed the elements by species names, so that we know which species the element was found in. For example:
perl ./step2/add_spp_name.pl s1.strained s1 s1.renamed.strained

cat *.renamed.strained > combined.renamed.strained
mmseqs easy-cluster combined.renamed.strained combined.renamed.strained.easyCluster tmp --min-seq-id 0.8 -c 0.8 --cov-mode 1 -v 2 --threads 10 --cluster-reassign

# Lastly, we annotated all the genomes using this library
earlGreyAnnotationOnly -g kaburagiaRhusicola.fa -s kaburagiaRhusicola -o output -t 20 -l combined.renamed.strained.easyCluster_rep_seq.fasta

```

## 3 - X chromosome assignment
Aphids have the XO sex determination system, where females have two copies of the X chromosomes, whereas males have one copy of the X. Therefore, we assigned X chromosomes in the genome assembly by estimating the sequencing depth of each chromosome using whole-genome sequencing data from males and females of *K. r. ovogallis*. As the X chromosomes were conserved in aphids (Li et al. 2020; Mathers et al. 2021; Li et al. 2023), we assigned the X chromosomes of other sumac gall aphid species based on the genome synteny with *K. r. ovogallis*. 

### 3.1 - Assigning X chromosomes of *K. r. ovogallis* based on the sequencing depth of males and females
We aligned reads from male and female *K. r. ovogallis* individuals to the reference genome. Chromosomes exhibiting half the sequencing depth in males compared with females were identified as X chromosomes. We then numbered the chromosomes in descending order of length. All sumac gall aphids have 10 autosomes (1 to 10) and 3 Xs (X1 to X3). 

<details>
<summary> Click to show the code for mapping </summary>
	
```
~/bin/bwa_v0.7.19/bwa index genome.fa
~/bin/bwa_v0.7.19/bwa mem -t 20 -o genome.sam ../../references/s5h1_plus_allBcomps.fa pe.1.fq.gz pe.2.fq.gz

samtools view -h -b -S genome.sam -o genome.bam --threads 20
samtools sort genome.bam -o genome.sorted.bam --threads 20
samtools index genome.sorted.bam
rm genome.sam genome.bam

# estimating the sequencing depth using mosdepth version 0.3.8
## first, we made sliding windows of 100kb and 50kb for each step
perl ./step2/length_distribution_V1.5.pl genome.fa
cut -f 1,2 genome.fa.length > genomeFile.txt
~/bin/bedtools_v2.31.0 makewindows -g genomeFile.txt -w 100000 -s 50000 > windows.bed

## then we estimated the sequencing depth of each window
mosdepth -b windows.bed -f genome.fa -n -t 20 --mapq 0 window_mapq0 genome.sorted.bam
mosdepth -b windows.bed -f genome.fa -n -t 20 --mapq 1 window_mapq1 genome.sorted.bam

# we also calculated the mapping rates
conda activate ~/miniconda3/envs/samtools_v1.15.1
## unmapped
echo "unmapped"
samtools view -c -f 4 --threads 10 genome.sorted.bam

## mapped & primary
echo "mapped"
samtools view -c -F 0x904 --threads 10 genome.sorted.bam
```

</details>

We plotted the distribution of sequencing depth in RStudio

<details>
<summary> <b><u>Click</b></u> to show the R code </summary>
	
```
depth_mapq0 <- read.table('window_mapq0.regions.bed.gz', header = F)
depth_mapq1 <- read.table('window_mapq1.regions.bed.gz', header = F)

depth_mapq0$V5 = "mapq0"
depth_mapq1$V5 = "mapq1"

depth <- rbind(depth_mapq0, depth_mapq1)

idx <- depth$V1 %in% names(key)
depth$V1[idx] <- key[depth$V1[idx]]

head(depth)
tail(depth)

#### boxplot - pick positive Bs and all As ####
depth_AB <- depth %>%
  mutate(
    chr = case_when(
      grepl("X", V1)      ~ "X",
      grepl("B$", V1)     ~ V1,
      TRUE                ~ "A"
    )
  )

depth_AB <- subset(depth_AB, chr %in% c("A", "X",
                                           "B1_s5h1_14B", "B1_s4_14B", 
                                           "B2_s5h1_15B", "B2_s5h1_16B", "B2_s5h1_17B",
                                           "B2_s5h1_18B"))

depth_AB$chr = factor(depth_AB$chr, levels = c("A", "X",
                                           "B1_s5h1_14B", "B1_s4_14B", 
                                           "B2_s5h1_15B", "B2_s5h1_16B", "B2_s5h1_17B",
                                           "B2_s5h1_18B"))


ggplot(depth_AB, aes(x = chr, y = V4, fill = V5)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.shape = NA
  ) + 
  coord_cartesian(ylim = c(0, 200)) + 
  scale_y_continuous(breaks = seq(0, 200, 50)) + # , expand = c(0, 0)
  theme_bw() +
  xlab("") +
  ylab("Depth") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.title = element_blank(),
    #legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.text.x  = element_blank(),
    #axis.ticks.x = element_blank()
  )
```

</details>

### 3.2 - Assigning X chromosomes of other species based on the genome synteny with *K. r. ovogallis*
We aligned assemblies to the genome of *K. r. ovogallis* using minimap2 version 2.28 and visualized using Asynt (https://github.com/simonhmartin/asynt). 

```
~/bin/minimap2-2.28_x64-linux/minimap2 -x asm20 -t 10 q.fa ref.fa | gzip > q_ref.paf.gz

samtools faidx q.fa
samtools faidx ref.fa


# we combined the two fai files and made a chr file so we could keep the order of chromosomes in the visualization step
cat q.fa.fai ref.fa.fai > q_ref.fa.fai
cat q_ref.fa.fai | perl -ne '@a=split; print "$a[0]\t$a[0]\t+\n";  '> q_ref_chr.txt 

```

In RStudio, we visualized the syntenies:
```
source("~/bin/asynt-master/asynt.R")

# load the output from minimap2
alignments <- import.paf("q_ref.paf.gz")

ref_data <- import.genome(fai_file="q_ref.fa.fai", chrom_file = "q_ref_chr.txt")
query_data <- import.genome(fai_file="q_ref.fa.fai", chrom_file = "q_ref_chr.txt")

# Remove short alignments
alignments_10kb <- subset(alignments, Rlen >= 10000 & Qlen >= 10000)

pdf("q_ref.pdf", width = 10, height = 10, useDingbats = FALSE) # Open a new pdf file
par(mar=c(1.5,1.5,0,0), xpd=NA)
plot.alignments.diagonal(alignments_10kb, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, reference_ori=ref_data$seq_ori, query_ori = query_data$seq_ori, no_reverse=FALSE, no_reorder=TRUE)
dev.off() # Close the file

```

The resulting dotplot looks like this:
<img width="1780" height="1770" alt="image" src="https://github.com/user-attachments/assets/dd2775f6-1916-45c8-8782-6349a8fd67d5" />


## 4 - Phylogeny reconstruction
### 4.1 - Ortholog inference
To understand the phylogenetic relationships of sumac gall aphids, we reconstructed the phylogeny using our species and 34 published aphid genomes (Table S2), including 1 species of Phylloxeridae (Family), 3 species of Adelgidae (Family), and 30 species of Aphididae (Family). We also used 11 Hemiptera species as outgroups, covering 3 species of Heteroptera (Suborder), 3 species of Auchenorrhyncha (Suborder), and 5 species of Sternorrhyncha (Suborder).

We first picked the longest isoform of each gene, then used OrthoFinder version 3.1.0 to identify orthologous groups. The amino acid sequences used can be found on Zenodo TBD.
```
orthofinder -f ./ -t 40 -a 5
```
### 4.2 - Phylogenetic tree construction
We used MAFFT version 7.526 with the L-INS-i algorithm to align amino acid sequences and Gblocks version 0.91b to remove regions that were not well aligned. Lastly, we concatenated the resulting alignments and constructed phylogenetic trees using IQ-Tree version 3.0.1 

```
mafft-linsi gene.fa 1> gene.fa.align 2> gene.fa.align.log
Gblocks gene.fa.align -t=p -b4=10 -b5=n -e=-gb
perl gblock_convert.pl gene.fa.align-gb > gene.GBlocks.fas

# iqtree
iqtree3 -s Orthogroups.GeneCount.tsv.0.80.concate.fasta -B 1000 -T 10 -o Beta,Dici,Pave,Trru,Rhpr,Onfa,Nilu,Sofu,Last,Erpe,Phso -redo
```

### 4.3 - Divergence time estimation
To estimate the divergence time, we used r8s version 1.81 with the penalized likelihood method and the TN algorithm based on fossil calibration points from Mathers et al. 2021 Molecular Biology and Evolution. 

In our case: 
```
mrca Reduviidae Trru Rhpr Onfa;
mrca AucHet Trin Rhpr Onfa Nilu Sofu Last;
mrca Sternorryncha Beta Dici Pave Erpe Phso Acpi AcJI Adab Adts Adco Apfa Apgo ApJI Brbr Chov Cice Drpl Davi Erla Hoco Maal Mecr Mesa Mype MyG6 Nefo Peni Peac Phhu Psba Rhma Rhpa Sehe Siav Simi Slle Teni Thtr Tusa s1 s10h1 s10h2 s11h1 s6h1 s6h2 s11h2 s3h1 s3h2 s4 s5h2 s5h1;
mrca Hemiptera Beta Dici Pave Erpe Phso Acpi AcJI Adab Adts Adco Apfa Apgo ApJI Brbr Chov Cice Drpl Davi Erla Hoco Maal Mecr Mesa Mype MyG6 Nefo Peni Peac Phhu Psba Rhma Rhpa Sehe Siav Simi Slle Teni Thtr Tusa s1 s10h1 s10h2 s11h1 s6h1 s6h2 s11h2 s3h1 s3h2 s4 s5h2 s5h1 Trin Rhpr Onfa Nilu Sofu Last;
fixage taxon=Reduviidae age=250;
fixage taxon=AucHet age=362;
fixage taxon=Sternorryncha age=356;
fixage taxon=Hemiptera age=386;
```
To run r8s, we used: 
```
r8s -b -f ./step4/r8s_ctl_file.modified.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > aphid_tree.txt
```

## 5 - Chromosomal rearrangement rate estimation
To quantitatively estimate chromosomal rearrangement rates, we calculated the amount of chromosomal rearrangements normalized by divergence time (per million years). 

### 5.1 - Divergence time between species pairs
To calculate rearrangement rates between species pairs, we first estimated divergence times between pairs of species. We used R to extract divergence times from the time tree from r8s. 

<details>
<summary> Click to show the R script to extract divergence time </summary>
	
```
# ref: https://rdrr.io/cran/ape/man/cophenetic.phylo.html
# ref: https://blog.phytools.org/2013/10/finding-edge-lengths-of-all-terminal.html

library(ape)
library(phytools)
library(reshape)
library(ggplot2)
library(gridExtra)
library(ggtree)
library(ggrepel)

tree <- ape::read.tree("./step5/aphid_tree.txt")
plotTree(tree)
edgelabels(tree$edge.length)

# pairwise distance
dist_m <- cophenetic.phylo(tree)

dist_c <- melt(dist_m, varnames = c("sp1", "sp2"), value.name = "value")
head(dist_c)

# the divergence time between the two species pairs equals to their distance divided by 2
dist_c$time = dist_c$value/2 
head(dist_c)

write.table(dist_c, file = "pairwise_time.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```

</details>

### 5.2 - Estimating intra-chromosomal rearrangement rates
To calculate the intra-chromosomal rearrangement rates, we first ran MCScanX version 1.0.0 and calculated the 

```
# step 0 - McScanX
## extract the 1-to-1 orthologs
perl ./step5/step0.pick_single_copy_genes.v20250915.pl Orthogroups.tsv s5h1 s1 combined.homology
# make gff file
cat s5h1.gff s1.gff > combined.gff
# run MCScanX
~/bin/MCScanX/MCScanX_h combined

# rates
# step1 - from .collinear to .bed file, labeling the start and end for each collinear block
perl ./step5/step1.MCScanX_collinearity_2_block_bed.pl combined.collinearity combined.gff dt sc combined.collinearity.bed

# step2 - orient the blocks from beginning to end based on sp1
perl ./step5/step2.orient_from_start2end.v20251030.pl combined.collinearity.bed combined.collinearity.bed.orient

# step3 - classify SVs
python step3.classify_rearrangements.v20251031.py -i combined.collinearity.bed.orient -o combined.collinearity.bed.orient.sv

```

### 5.3 - Estimating inter-chromosomal rearrangement rates

```
# inter-chromosomal rearrangement rate estimation - v20260411
# step5 - calculate inter-chr changes
perl ./step5/step5.inter_chr.v20260414.pl combined.collinearity.bed.orient.sv_classified.tsv s5h1.chr.list s1.chr.list step5.inter_v20260411

# step6 - calculate rearrangement density
perl ./step5/step6.fussion_fission.v20260407.pl step5.inter_v20260411 ./step5/X.list step6.output_v20260411
```

### 5.4 - Visualization of rates

In R, we used the following commands to visualize the intra-chromosomal rearrangement rates:

<details>
<summary> Click to show the R script for visualization and statistical tests </summary>

```
# 1 - intra-chr rates
# 1.1 - read tables
library(ggplot2)
library(ggrepel)

intra_Adelgidae <- read.table('./step5/Adelgidae.intra_sv.list.rates', header = T)
intra_Adelgidae$group = "Adelgidae"

intra_Aphidinae <- read.table('./step5/Aphidinae.intra_sv.list.rates', header = T)
intra_Aphidinae$group = "Aphidinae"

intra_Aphidini <- read.table('./step5/Aphidini.intra_sv.list.rates', header = T)
intra_Aphidini$group = "Aphidini"

intra_Fordini <- read.table('./step5/Fordini.intra_sv.list.rates', header = T)
intra_Fordini$group = "Fordini"

intra_Macrosiphini <- read.table('./step5/Macrosiphini.intra_sv.list.rates', header = T)
intra_Macrosiphini$group = "Macrosiphini"

# outgroups
intra_Delphacinae <- read.table('./step5/Delphacinae.intra_sv.list.rates', header = T)
intra_Delphacinae$group = "Delphacinae"

intra_Psylloidea <- read.table('./step5/Psylloidea.intra_sv.list.rates', header = T)
intra_Psylloidea$group = "Psylloidea"

intra_Triatominae <- read.table('./step5/Triatominae.intra_sv.list.rates', header = T)
intra_Triatominae$group = "Triatominae"

# exclude a few species
intra_combined_ingroup <- rbind(intra_Adelgidae,intra_Aphidini, intra_Fordini, intra_Macrosiphini)

intra_noSimi <- intra_combined_ingroup[!grepl("Simi", intra_combined_ingroup$spp), ]

intra_noSimi_noChov <- intra_noSimi[!grepl("Chov", intra_noSimi$spp), ]
intra_noSimi_noChov_nos4 <- intra_noSimi_noChov[!grepl("s4", intra_noSimi_noChov$spp), ]

# divergence time between 10 to 35 MYA
intra_noSimi_noChov_nos4_mid <- subset(intra_noSimi_noChov_nos4, time >= 10 & time <= 35)

# 1.2 estimate statistics
library(dplyr)

# the rate is number of intra-chr rearrangements per Mb per MY
intra_noSimi_noChov_nos4 %>% 
  group_by(group) %>%
  summarise(
    median = median(rates, na.rm = TRUE),
    Q1 = quantile(rates, 0.25, na.rm = TRUE), 
    Q3 = quantile(rates, 0.75, na.rm = TRUE)
  )

# 1.3 - violin plot
intra_noSimi_noChov_nos4_mid$group = factor(intra_noSimi_noChov_nos4_mid$group, levels = c("Adelgidae", "Fordini", "Aphidini", "Macrosiphini"))
p_intra <- ggplot(data = intra_noSimi_noChov_nos4_mid, aes(x = group, y = rates, group = group, color = group)) +
			geom_violin(aes(color = group)) +
			geom_point(position = position_jitter(seed = 1, width = 0.2, height = 0)) +
			theme(legend.position = "none") +
			theme_bw() +
			geom_label_repel(aes(label = spp),
				max.overlaps = 200,
				box.padding   = 0.35,
				point.padding = 0.5,
				segment.color = 'grey50') +
			theme(
				plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
				legend.position= "none",
				strip.background=element_blank(),
				axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
				axis.title.x = element_text(size = 15), 
				axis.text.y = element_text(size = 15), 
				axis.title.y = element_text(size = 15)
			) +
			xlab('Taxonomic groups') + ylab('Intra-chr rearrangement rates') # + ggtitle('Sequencing depth')

# 1.4 - kruskal-wallis tests
kruskal.test(rates ~ group, data = intra_noSimi_noChov_nos4_mid)

library(FSA)

dunnTest(
  rates ~ group,
  data = intra_noSimi_noChov_nos4_mid,
  method = "bh"
)

# 2 - inter-chr rates
library(ggplot2)
library(ggrepel)

inter_Adelgidae <- read.table('./step5/Adelgidae.inter_sv.list.rates', header = T)
inter_Adelgidae$group = "Adelgidae"
head(inter_Adelgidae)

inter_Aphidinae <- read.table('./step5/Aphidinae.inter_sv.list.rates', header = T)
inter_Aphidinae$group = "Aphidinae"
head(inter_Aphidinae)

inter_Aphidini <- read.table('./step5/Aphidini.inter_sv.list.rates', header = T)
inter_Aphidini$group = "Aphidini"
head(inter_Aphidini)

inter_Fordini <- read.table('./step5/Fordini.inter_sv.list.rates', header = T)
inter_Fordini$group = "Fordini"
head(inter_Fordini)

inter_Macrosiphini <- read.table('./step5/Macrosiphini.inter_sv.list.rates', header = T)
inter_Macrosiphini$group = "Macrosiphini"
head(inter_Macrosiphini)

# outgroups
inter_Delphacinae <- read.table('./step5/Delphacinae.inter_sv.list.rates', header = T)
inter_Delphacinae$group = "Delphacinae"
head(inter_Delphacinae)


inter_Psylloidea <- read.table('./step5/Psylloidea.inter_sv.list.rates', header = T)
inter_Psylloidea$group = "Psylloidea"
head(inter_Psylloidea)

inter_Triatominae <- read.table('./step5/Triatominae.inter_sv.list.rates', header = T)
inter_Triatominae$group = "Triatominae"
head(inter_Triatominae)

inter_combined <- rbind(inter_Adelgidae,inter_Aphidini, inter_Fordini, inter_Macrosiphini, 
                        inter_Delphacinae, inter_Psylloidea, inter_Triatominae)

# exclude a few species
inter_noSimi <- inter_combined[!grepl("Simi", inter_combined$spp), ]
inter_noSimi_noChov <- inter_noSimi[!grepl("Chov", inter_noSimi$spp), ]
inter_noSimi_noChov_nos4 <- inter_noSimi_noChov[!grepl("s4", inter_noSimi_noChov$spp), ]

# divergence time between 10 - 35 MYA
inter_noSimi_noChov_nos4_mid <- subset(inter_noSimi_noChov_nos4, time >= 10 & time <= 35)

# 2.2 - estimate statistics
# the inter-chr rate is the number of intra-chr rearrangements per Mb per MY
inter_noSimi_noChov_nos4 %>% 
  group_by(group) %>%
  summarise(
    median = median(inter_rate_perMY, na.rm = TRUE),
    Q1 = quantile(inter_rate_perMY, 0.25, na.rm = TRUE), 
    Q3 = quantile(inter_rate_perMY, 0.75, na.rm = TRUE)
  )

# 2.3 - violin plot
inter_noSimi_noChov_nos4_mid$group = factor(inter_noSimi_noChov_nos4_mid$group, 
                                            levels = c("Delphacinae", "Adelgidae", 
                                                       "Fordini", "Aphidini", "Macrosiphini"))

(
  p_inter <- ggplot(data = inter_noSimi_noChov_nos4_mid, aes(x = group, y = inter_rate_perMY, group = group, color = group)) +
  geom_violin(aes(color = group)) + 
  geom_point(position = position_jitter(seed = 1, width = 0.2, height = 0)) +
  theme(legend.position = "none") + 
  theme_bw() + 
        # geom_label_repel(aes(label = spp),
        #            max.overlaps = 200,
        #            box.padding   = 0.35,
        #            point.padding = 0.5,
        #            segment.color = 'grey50') +
  theme(#strip.text.y=element_text(angle=180, size=15),
        #strip.text.x=element_text(angle=180, size=8),
        #legend.position=c(.99,.01), legend.justification=c(1,0),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
        # plot.title = element_text(size = 40, )
        legend.position= "none",
        strip.background=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15, color = "black"),
        axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        #axis.text.y=element_blank(), axis.ticks.y=element_blank()
        ) + 
    xlab('Taxonomic groups') + ylab('Inter-chromosomal rearrangement rates per MY') # + ggtitle('Sequencing depth')
)

# 2.4 - kruskal-wallis test
library(FSA)

kruskal.test(inter_rate_perMY ~ group, data = inter_noSimi_noChov_nos4)

dunnTest(
  inter_rate_perMY ~ group,
  data = inter_noSimi_noChov_nos4,
  method = "bh"
)
```

</details>


## 6 - Improving the B chromosome assembly
To improve the assembly of B chromosomes, we joined B scaffolds based on four lines of evidence: (1) hifiasm assemblies of B-derived reads, (2) assembly graphs of Verkko, (3) HiC contact maps, and (4) genome syntenies. A joining was made if three out of four lines of evidence supported such connection. 

### 6.1 - Joining based on hifiasm assemblies of B-derived reads
We aligned HiFi reads to the original assembly and extracted HiFi reads that did not map to autosomes or Xs. We assigned these reads as putative B-derived reads. We only assembled these B-derived reads using hifiasm to simplify the assembly process. We aligned the original B scaffolds to the B-derived scaffolds. If the B-derived assembly connected two B scaffolds from the original assembly, we treated this connection as supported by B-derived hifiasm evidence. 

```
minimap2 -a genome.fa -o genome.sam -t 20 -x map-hifi hifi.filt.fastq.gz
samtools view -h -b -S genome.sam -o genome.bam --threads 20
samtools sort genome.bam -o genome.sorted.bam --threads 20
samtools index genome.sorted.bam

# example: h1tg000021l ... h2tg000129l are the scaffolds that are not scaffold-level
samtools view -bS genome.sorted.bam -o B_reads.bam h1tg000021l h1tg000016l h1tg000020l h1tg000046l h1tg000011l h1tg000025l h1tg000018l h1tg000026l h1tg000034l h1tg000031l h1tg000022l h1tg000029l h1tg000028l h1tg000040l h1tg000038l h1tg000043l h1tg000037l h1tg000041l h1tg000017l h1tg000035l h1tg000023l h1tg000062l h1tg000036l h1tg000058l h1tg000068l h1tg000047l h1tg000056l h1tg000027l h1tg000061l h1tg000033l h1tg000039l h1tg000059l h1tg000052l h1tg000067l h1tg000072l h1tg000070l h1tg000014l h1tg000054l h1tg000051l h1tg000032l h1tg000045l h1tg000071l h1tg000079l h1tg000055l h1tg000064l h1tg000065l h1tg000049l h1tg000084l h1tg000030l h1tg000085l h1tg000050l h1tg000075l h1tg000063l h1tg000078l h1tg000042l h1tg000073l h1tg000066l h1tg000074l h1tg000091l h1tg000069l h1tg000081l h1tg000044l h1tg000060l h1tg000057l h1tg000076l h1tg000053l h1tg000088l h1tg000024l h1tg000089l h1tg000082l h1tg000092l h1tg000090l h1tg000087l h1tg000083l h1tg000048l h1tg000086l h1tg000077l h1tg000080l h2tg000001l h2tg000057l h2tg000012l h2tg000003l h2tg000009l h2tg000026l h2tg000039l h2tg000029l h2tg000014l h2tg000033l h2tg000031l h2tg000046l h2tg000059l h2tg000069l h2tg000061l h2tg000017l h2tg000025l h2tg000086l h2tg000011l h2tg000064l h2tg000018l h2tg000038l h2tg000034l h2tg000036l h2tg000068l h2tg000021l h2tg000053l h2tg000037l h2tg000035l h2tg000098l h2tg000058l h2tg000041l h2tg000056l h2tg000051l h2tg000062l h2tg000072l h2tg000042l h2tg000043l h2tg000075l h2tg000063l h2tg000048l h2tg000097l h2tg000084l h2tg000013l h2tg000079l h2tg000032l h2tg000047l h2tg000085l h2tg000082l h2tg000078l h2tg000045l h2tg000066l h2tg000090l h2tg000077l h2tg000054l h2tg000016l h2tg000091l h2tg000081l h2tg000067l h2tg000092l h2tg000055l h2tg000107l h2tg000112l h2tg000071l h2tg000116l h2tg000060l h2tg000106l h2tg000094l h2tg000052l h2tg000050l h2tg000065l h2tg000099l h2tg000040l h2tg000088l h2tg000101l h2tg000117l h2tg000074l h2tg000105l h2tg000093l h2tg000073l h2tg000108l h2tg000118l h2tg000110l h2tg000124l h2tg000049l h2tg000125l h2tg000076l h2tg000119l h2tg000104l h2tg000120l h2tg000083l h2tg000089l h2tg000115l h2tg000044l h2tg000070l h2tg000114l h2tg000095l h2tg000096l h2tg000087l h2tg000080l h2tg000113l h2tg000102l h2tg000109l h2tg000103l h2tg000019l h2tg000100l h2tg000123l h2tg000030l h2tg000126l h2tg000130l h2tg000128l h2tg000121l h2tg000122l h2tg000127l h2tg000131l h2tg000111l h2tg000129l --threads 20

samtools fastq B_reads.bam | gzip > B_reads.fastq.gz

```

Then assemble the reads and find joins that connect two B-derived scaffolds using minimap2:
```
~/bin/minimap2-2.28_x64-linux/minimap2 -o ccs --dual-scaf --telo-m CCTAA --primary -t 40 B_reads.fastq.gz
awk '/^S/{print ">"$2;print $3}' ccs.p_ctg.gfa > ccs.p_ctg.fa

~/bin/minimap2-2.28_x64-linux/minimap2 -c -x asm5 -t 20 genome.B.fa assembly.fasta > verkko_vs_hifiasm.asm5_CIGAR.paf
```

### 6.2 - Joining based on the assembly graphs of Verkko
We assembled B-derived reads with Verkko version 2.2.1. We then aligned the original hifiasm B scaffolds to the assembly graph using GraphAligner version 1.0.20. If the assembly graph connected two B scaffolds, we treated the connection as supported by the Verkko assembly graph

```
verkko -d run1_asm --hifi B_reads.fastq.gz --telomere-motif CCTAA
GraphAligner -g assembly.homopolymer-compressed.gfa -f genome.B.fa -a genome.B.graph_aligner.dbg.gaf -x dbg --threads 10
cat genome.B.graph_aligner.dbg.gaf | awk '($11 > 100000)&&($12 >= 60)' | sort -k3,3n > genome.B.graph_aligner.dbg.gaf.filtered

```

### 6.3 - Joining based on HiC contact maps
Thirdly, if HiC maps suggested contacts between two scaffolds, we treated the joining as supported by HiC signals. The steps we used to make HiC contact maps can be found in the Manual curation section [2.4](#24---manual-curation). 

### 6.4 - Joining based on genome syntenies
Lastly, if other complete Bs showed synteny of two fragmented B scaffolds, we treated the joining as supported by synteny. The steps we used to make genome synteny can be found in section [3.2](#32---assigning-x-chromosomes-of-other-species-based-on-the-genome-synteny-with-k-r-ovogallis). 

## 7 - *dN/dS* and gene expression of B-linked genes
### 7.1 - Measuring *dN/dS*
To measure the amount of non-synonymous substitutions versus synonymous substitutions, we first extracted the single-copy orthologs among the five species: *K. r. ovogallis*, *K. r. ovatirhusicola*, *K. r. rhusicola*, *K. r. ensigallis*, *K. r. sp.*, and *Floraphis meitanensis*. To measure the differences between A-linked genes and B-linked genes, we treated each B chromosome as an independent sample. So each gene must be single-copy on either A chromosomes or a B chromosome. Then we extracted the CDS of these single-copy orthologs, performed codon-based alignments (https://github.com/lyy005/codon_alignment), and calculated *dN/dS* ratios. We used *Floraphis meitanensis* as the outgroup. The *dN/dS* of each gene copy was measured in comparison to the gene copy of *Floraphis meitanensis* using KaKs_Calculator version 3.0 (https://github.com/Chenglin20170390/KaKs_Calculator-3.0). 

```
# first, ran codon alignments
perl ~/bin/codon_alignment_v1.5/codon_alignment.pl gene.CDS.fa 1 0

# Then, convert alignments to .axt format
perl ./step7/aln2axt.pl gene.CDS.fa.nt.cleanup.aln

# lastly, calculated dN/dS ratios using KaKs calculator
~/bin/KaKs_Calculator_v3.0/bin/KaKs -i gene.CDS.fa.nt.cleanup.aln.pw.axt -o gene.kaks
```

### 7.2 - Estimating TPM 
To estimate the expression levels of genes on As (autosomes and Xs) and Bs, we measured the Transcripts Per Million (TPM) for each gene using Salmon version 2.3.2. 

We first extracted the longest CDS of each gene:
```
perl ./step2/find_longest_protein_Braker_v20210403b.braker.pl braker.cds.fa braker.cds.fa.longest

```
Then we used salmon to measure TPM values. 
```
salmon index --keepDuplicates -t braker.cds.fa.longest -i salmon_index_keepDup -p 10
salmon quant -i salmon_index_keepDup -l A -1 pe.1.fq.gz -2 pe.2.fq.gz -p 10 -o salmon_quant_keepDup
```

## 8 - 3D organization of B chromosomes
### 8.1 - TAD
To explore the 3D organization of Bs, we followed the steps of HiCExplorer version 3.7.6. 

First, we mapped HiC reads to aphid genomes with 1 copy of each B type. To avoid multi-mapped reads, we removed reads with MAPQ < 15. Then we corrected the HiC matrix, merged to 30kb resolution, and identified TADs:
```
~/bin/bwa_v0.7.19/bwa index genome.fa
~/bin/bwa_v0.7.19/bwa mem -A 1 -B 4 -E 50 -L 0 -t 10 genome.fa mate_R1.fastq.gz 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam
~/bin/bwa_v0.7.19/bwa mem -A 1 -B 4 -E 50 -L 0 -t 10 genome.fa mate_R2.fastq.gz 2>>mate_R2.log | samtools view -Shb - > mate_R2.bam

hicBuildMatrix --samFiles mate_R1.bam mate_R2.bam --binSize 10000 --restrictionSequence GATC --danglingSequence GATC --restrictionCutFile cut_sites.bed --threads 20 --inputBufferSize 100000 --outBam ./hicMatrix/hic.bam -o ./hicMatrix/hic_matrix.h5 --QCfolder ./hicMatrix/hicQC --minMappingQuality 15 > step2.1.log 2> step2.2.log &

hicCorrectMatrix diagnostic_plot -m hic_matrix.h5 -o hic_matrix.diagnostic_plot.pdf
hicCorrectMatrix correct -m hic_matrix.h5 -o hic_corrected.h5 --correctionMethod KR

hicMergeMatrixBins --matrix hic_corrected.h5 --outFileName hic_corrected_30kb.h5 --numBins 3
hicFindTADs -m hic_corrected_30kb.h5 --outPrefix hic_corrected_30kb.findTADs --numberOfProcessors 16 --correctForMultipleTesting fdr

# Plot TADs
hicPlotTADs --tracks ./step8/hic_tracks.corrected.ini --region s5h1_14B:0-13788911 -o hic_corrected_30kb.hicPlotTADs.s5h1_14B.pdf
hicPlotTADs --tracks ./step8/hic_tracks.corrected.ini --region s5h1_15B:0-12963132 -o hic_corrected_30kb.hicPlotTADs.s5h1_15B.pdf

```
### 8.2 - trans-interactions
We calculated the average of interactions between Bs and As (autosomes and Xs) using cooltools version 0.7.1. 

```
hicConvertFormat --matrices hic_matrix.h5 --inputFormat h5 --outputFormat cool --outFileName hic_matrix.cool

conda activate cooltools_v0.7.1
cooler balance hic_matrix.cool
cooltools expected-trans hic_matrix.cool -p 10 -o trans_contacts.tsv
```

We visualized the trans-contacts in R:

<details>
<summary> Click to show the R script for visualization </summary>

```
df <- read.table('trans_contacts.tsv', header = T)

head(df)

# select only the chrs
selected_chr <- c(
  "s5h1_1", "s5h1_2", "s5h1_3", "s5h1_4", "s5h1_5", "s5h1_6", "s5h1_7", "s5h1_8", "s5h1_9", "s5h1_10", 
  "s5h1_X1", "s5h1_X2", "s5h1_X3", "s5h1_14B", "s5h1_15B"
)


df_subset <- df %>%
  filter(
    region1 %in% selected_chr,
    region2 %in% selected_chr
  ) %>%
  distinct(region1, region2, .keep_all = TRUE)

df_subset


# Create full symmetric pair matrix
mat <- bind_rows(
  df_subset %>% select(r1 = region1, r2 = region2, val = balanced.avg),
  df_subset %>% select(r1 = region2, r2 = region1, val = balanced.avg)
) %>%
  distinct(r1, r2, .keep_all = TRUE) %>%
  pivot_wider(names_from = r2, values_from = val) %>%
  column_to_rownames("r1") %>%
  as.matrix()

# Define your exact desired order
target_order <- c(
  "s5h1_1", "s5h1_2", "s5h1_3", "s5h1_4", "s5h1_5", 
  "s5h1_6", "s5h1_7", "s5h1_8", "s5h1_9", "s5h1_10", 
  "s5h1_X1", "s5h1_X2", "s5h1_X3", "s5h1_14B", "s5h1_15B"
)

# Subset to keep only regions present in your current matrix
valid_order <- target_order[target_order %in% rownames(mat)]

# Reorder matrix rows and columns
mat_ordered <- mat[valid_order, valid_order]

# Set the diagonal to NA
diag(mat_ordered) <- NA

pheatmap(
  mat_ordered,
  angle_col = "45",
  cluster_rows = F,
  cluster_cols = F,
  scale = "none",
  fontsize_row = 8,
  fontsize_col = 8, 
  na_col = "grey90",
)
```
</details>

## 9 - Methylation in B chromosomes
To identify DNA methylation, we used pbmm2 version 1.17.0 (https://github.com/PacificBiosciences/pbmm2) and aligned_bam_to_cpg_scores version 3.0 (https://github.com/PacificBiosciences/pb-CpG-tools). We combined both haplotypes of *K. r. ensigallis* so that the sequencing depth of autosomes and Xs were similar to that of Bs. 

### 9.1 - Methylation detection
```
pbmm2 index genome.fa genome.mmi
pbmm2 align genome.mmi hifi_reads.bam aligned.bam --preset CCS --sort --log-level DEBUG --num-threads 20

aligned_bam_to_cpg_scores --bam aligned.bam --output-prefix cpg_denovo_mapq20 --threads 20 --min-mapq 20 --modsites-mode denovo

# keep CpGs with coverage >= 5
less -S cpg_denovo_mapq20.combined.bed.gz | awk '$7 >= 5' > cpg_denovo_mapq20.coverage5.bed
```

### 9.2 - CpG sites in the genome
Then we measured sequencing depth from the output BAM file, and only kept sites with sequencing depth >= 5. 

```
mosdepth -t 20 window_mapq20 --mapq 20 aligned.bam

# mosdepth cov >= 5 bed
less -S window_mapq20.per-base.bed.gz | awk '$4 >= 5' > window_mapq20.per-base.coverage5.bed

# CpGs in the genome were extracted
perl ./step9/make_CpG_BEDs.pl genome.fa genome.bed

# the output is the BED with CpGs in genomes and with coverage >= 5 (and mapq >= 20)
~/bin/bedtools_v2.31.0 intersect -a genome.bed -b window_mapq20.per-base.coverage5.bed -u > genome_CpG_AND_coverage5.bed

```

### 9.3 - Identifying exon, intron, and intergenic regions
We combined the GFF files of both hap1 and hap2 and extracted exon, intron, and intergenic regions. Then we used bedtools to extract the CpGs in each region.

<details>
<summary> Click to show the code for extracting CpGs in different regions </summary>

```
cat hap1.gff hap2.gff > genome.gff

awk '$3 == "exon"' genome.gff | perl -ne '@a=split; $a[3] = $a[3] -1; $b=join "\t", @a; print "$a[0]\t$a[3]\t$a[4]\n"; ' > annotation.exon.bed

# find gene BEDs
awk '$3 == "gene"' genome.gff | perl -ne '@a=split; $a[3] = $a[3] -1; $b=join "\t", @a; print "$a[0]\t$a[3]\t$a[4]\n"; ' > annotation.gene.bed

# merge windows
## exon
~/bin/bedtools_v2.31.0 sort -i annotation.exon.bed > annotation.exon.bed.sorted
~/bin/bedtools_v2.31.0 merge -i annotation.exon.bed.sorted > annotation.exon.bed.sorted.merged

# gene
~/bin/bedtools_v2.31.0 sort -i annotation.gene.bed > annotation.gene.bed.sorted
~/bin/bedtools_v2.31.0 merge -i annotation.gene.bed.sorted > annotation.gene.bed.sorted.merged

# HiFi CpG
## exons
~/bin/bedtools_v2.31.0 intersect -a cpg_denovo_mapq20.coverage5.bed -b annotation.exon.bed.sorted.merged -u > hifi_CpG.exon.bed

## gene
~/binbedtools_v2.31.0 intersect -a cpg_denovo_mapq20.coverage5.bed -b annotation.gene.bed.sorted.merged -u > hifi_CpG.gene.bed

## intron: within gene but not in exon
~/bin/bedtools_v2.31.0 intersect -a hifi_CpG.gene.bed -b annotation.exon.bed.sorted.merged -v > hifi_CpG.intron.bed

## intergenic: not in gene
~/bin/bedtools_v2.31.0 intersect -a cpg_denovo_mapq20.coverage5.bed -b annotation.gene.bed.sorted.merged -v > hifi_CpG.intergenic.bed

# genome CpGs
## exons
~/bin/bedtools_v2.31.0 intersect -a genome_CpG_AND_coverage5.bed -b annotation.exon.bed.sorted.merged -u > genome_CpG.exon.bed

## gene
~/bin/bedtools_v2.31.0 intersect -a genome_CpG_AND_coverage5.bed -b annotation.gene.bed.sorted.merged -u > genome_CpG.gene.bed

## intron: within gene but not in exon
~/bin/bedtools_v2.31.0 intersect -a genome_CpG.gene.bed -b annotation.exon.bed.sorted.merged -v > genome_CpG.intron.bed

## intergenic
~/bin/bedtools_v2.31.0 intersect -a genome_CpG_AND_coverage5.bed -b annotation.gene.bed.sorted.merged -v > genome_CpG.intergenic.bed

for i in exon gene intron intergenic
do
 echo $i
 cut -f 1 genome_CpG.$i\.bed | sort | uniq -c > genome_CpG.$i\.each_chr.table
done
```

</details>

### 9.4 - Summarization and visualization of CpGs

<details>
<summary> Click to show the R script </summary>

```
# load libraries
library(data.table)
library(ggplot2)

# define types of chromosomes
A_chr <-  c("s11h1_1", "s11h1_2", "s11h1_3", "s11h1_4", "s11h1_5", "s11h1_6", "s11h1_7", "s11h1_8", "s11h1_9", "s11h1_10", 
            "s11h2_1", "s11h2_2", "s11h2_3", "s11h2_4", "s11h2_5", "s11h2_6", "s11h2_7", "s11h2_8", "s11h2_9", "s11h2_10")

X_chr <- c("s11h1_X1", "s11h1_X2", "s11h1_X3", "s11h2_X1", "s11h2_X2", "s11h2_X3")

B_chr <- c("s11h1_14B", "s11h1_15B", 
  "s11h1_16B", "s11h1_17B", "s11h1_18B", 
  "s11h1_19B", "s11h1_20B", "s11h1_21B", 
  "s11h1_31B", "s11h1_32B", "s11h1_33B", 
  "s11h1_34B", "s11h1_35B")

selected_chr <- c(A_chr, X_chr, B_chr)

# this is the code for introns, but we can make one for exons and intergenic regions
## 1 - number of CpGs detected by HiFi reads
meth_cov5_intron <- fread("hifi_CpG.intron.bed",
  header = F,
  select = c(1, 2, 4, 6), 
  col.names = c(
    "chrom",
    "begin",
    "mod_score",
    "cov"
  )
)

head(meth_cov5_intron)

# keep sites with on chrs
meth_cov5_intron_chr <- meth_cov5_intron[chrom %in% selected_chr]

meth_cov5_intron_chr[, group := fcase(
  chrom %in% A_chr, "A",
  chrom %in% X_chr, "X",
  chrom %in% B_chr, "B"
)]

meth_cov5_intron_chr$group = factor(meth_cov5_intron_chr$group, levels = c("A", "X", "B"))

## 2 - number of CpGs on the genome. This number is used later for calculating the percentage of methylated sites
intron_CpGs <- fread(
  "genome_CpG.intron.each_chr.table",
  select = 1:2,
  col.names = c("total_CpGs", "chr")
)

intron_CpGs <- intron_CpGs[chr %in% selected_chr]

# Assign chromosome groups
intron_CpGs[, group := fcase(
  chr %in% A_chr, "A",
  chr %in% X_chr, "X",
  chr %in% B_chr, "B"
)]

# Sum of CpGs for A, X, and B
group_intron_CpGs <- intron_CpGs[, .(
  total_CpGs_chr = sum(total_CpGs)
), by = group]

group_intron_CpGs


# 3 - plotting
```{r }
# Define bin boundaries
# breaks <- c(seq(0, 94, by = 5), 100)
breaks <- seq(0, 100, by = 10)

# Assign CpGs to bins
meth_cov5_intron_chr[, bin := cut(
  mod_score,
  breaks = breaks,
  include.lowest = TRUE,
  right = FALSE
)]

head(meth_cov5_intron_chr)

meth_cov5_intron_chr_hist <- meth_cov5_intron_chr[, .(
  N = .N
), by = .(group, bin)]

# Show the actual ranges and counts
print(meth_cov5_intron_chr_hist)


# add the total number of CpGs for A, X, B
meth_cov5_intron_chr_hist <- merge(
  meth_cov5_intron_chr_hist,
  group_intron_CpGs,
  by = "group"
)


head(meth_cov5_intron_chr_hist)


meth_cov5_intron_chr_hist$methylated_per_CpG = meth_cov5_intron_chr_hist$N / meth_cov5_intron_chr_hist$total_CpGs_chr

head(meth_cov5_intron_chr_hist)

meth_cov5_intron_chr_hist$group <- factor(meth_cov5_intron_chr_hist$group, levels = c("A", "X", "B"))

(
  p_meth_density_intron <-
    ggplot(meth_cov5_intron_chr_hist,
           aes(x = bin,
               y = methylated_per_CpG,
               color = group,
               group = group
               )
           ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    scale_color_manual(values = c(
      "A" = "#F8766D",
      "X" = "#619CFF",
      "B" = "#00BA38")
      ) +
    scale_y_log10(
      breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1),
      labels = scales::label_percent()
      ) +
    theme_bw() + 
    theme(# strip.text.y=element_text(angle=180, size=15),
        # strip.text.x=element_text(angle=180, size=8),
        #legend.position=c(.99,.01), legend.justification=c(1,0),
        # legend.position= "none",
        strip.background=element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.x = element_text(color = "black", size = 15, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 20, face = "bold", color = "black", hjust = 0.5)
        #plot.title = element_blank()
        #axis.text.y=element_blank(), axis.ticks.y=element_blank()
        ) + 
    xlab('mod_score interval') + ylab('Percentage of CpGs') + ggtitle('Intron')
)

```

</details>

CpGs on each type of chromosome: 

<details>
<summary> Click to show R script </summary>

```
# hifi plot
meth_cov5_exon_chr$region <- "exon"
meth_cov5_intron_chr$region <- "intron"
meth_cov5_intergenic_chr$region <- "intergenic"

hifi <- rbind(meth_cov5_exon_chr, meth_cov5_intron_chr, meth_cov5_intergenic_chr)

head(hifi)

# Count CpGs by source, chromosome group, and region
hifi_summary <- hifi[, .(N = .N), by = .(group, region)]

head(hifi_summary)

hifi_summary$group = factor(hifi_summary$group, levels = c("A", "X", "B"))

hifi_summary$region = factor(hifi_summary$region, levels = c("intergenic", "intron", "exon"))

(
  p_CpG_number <- ggplot(hifi_summary, aes(x = region, y = N / 1e6, fill = group)) +
  geom_col(color = "black", linewidth = 0.2) +
    scale_y_continuous(breaks = seq(0, 30, 5)) + 
  scale_fill_manual(values = c(
    "A" = "#F8766D",
    "X" = "#619CFF",
    "B" = "#00BA38"
  )) +
    theme_bw() + 
    theme(# strip.text.y=element_text(angle=180, size=15),
        # strip.text.x=element_text(angle=180, size=8),
        #legend.position=c(.99,.01), legend.justification=c(1,0),
        # legend.position= "none",
        strip.background=element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.x = element_text(color = "black", size = 15, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(color = "black", size = 15), 
        axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 20, face = "bold", color = "black", hjust = 0.5)
        #plot.title = element_blank()
        #axis.text.y=element_blank(), axis.ticks.y=element_blank()
        ) + 
    xlab('Type') + ylab('Number of CpGs (Mb)') + ggtitle('Number of CpGs')
)
```

</details>

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

### 2.2 - TRF
On top of telomeres found by seqtk, TRF also searches for telomeres in the middle or close to the end of each scaffold:
```
conda activate trf_v4.10.0
trf genome.fa 2 7 7 80 10 50 500 -f -d -m
```

## 3 - Characterize Bs in aphids
We also characterized other features of the genomes and the Bs



### 3.2 - Annotation of repeats
To annotate repetitive content, we used EarlGrey. 

We first annotated each genome separately:
```
conda activate earlgrey_v6.3.2
earlGrey -g kaburagiaRhusicola.fa -s kaburagiaRhusicola -r arthropoda -o output -t 20 -c yes -e yes -d yes
```

Then we combined all the .strained outputs:
```
for i in s10h1 s10h2 s11h1 s11h2 s1 s3h1 s3h2 s4 s5h1 s5h2 s6h1 s6h2
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





