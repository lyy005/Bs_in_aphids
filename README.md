# Bs_in_aphids
Pipeline for identifying supernumerary B chromosomes in sumac gall aphids


## 1 - Improving the assembly of Bs
### 1.1 - Extracting B-derived reads
To simplify the B chromsome assembly process, we first identified autosomes and Xs based on sequencing homology

Use SAMtools to keep reads that did not align to autosomes and Xs
```
samtools view -bS genome.sorted.bam -o B_reads.bam s5h1_14F s5h1_15F s5h1_16F s5h1_17F s5h1_18F s5h1_19F s5h1_20F s5h1_21F s5h1_22F s5h1_23F s5h1_24F s5h1_25F s5h1_26F s5h1_27F s5h1_28F s5h1_29F s5h1_30F s5h1_31F s5h1_32F s5h1_33F s5h1_34F s5h1_35F s5h1_36F s5h1_37F s5h1_38F s5h1_39F s5h1_40F s5h1_41F s5h1_42F s5h1_43F s5h1_44F s5h1_45F s5h1_46F s5h1_47F s5h1_48F s5h1_49F s5h1_50F s5h1_51F s5h1_52F s5h1_53F s5h1_54F s5h1_55F s5h1_56F s5h1_57F s5h1_58F s5h1_59F s5h1_60F s5h1_61F s5h1_62F s5h1_63F s5h1_64F s5h1_65F s5h1_66F s5h1_67F s5h1_68F s5h1_69F s5h1_70F s5h1_71F s5h1_72F s5h1_73F s5h1_74F s5h1_75F s5h1_76F s5h1_77F s5h1_78F s5h1_79F s5h1_80F s5h1_81F s5h1_82F s5h1_83F s5h2_14F s5h2_15F s5h2_16F s5h2_17F s5h2_18F s5h2_19F s5h2_20F s5h2_21F s5h2_22F s5h2_23F s5h2_24F s5h2_25F s5h2_26F s5h2_27F s5h2_28F s5h2_29F s5h2_30F s5h2_31F s5h2_32F s5h2_33F s5h2_34F s5h2_35F s5h2_36F s5h2_37F s5h2_38F s5h2_39F s5h2_40F s5h2_41F s5h2_42F s5h2_43F s5h2_44F s5h2_45F s5h2_46F s5h2_47F s5h2_48F s5h2_49F s5h2_50F s5h2_51F s5h2_52F s5h2_53F s5h2_54F s5h2_55F s5h2_56F --threads 20
```

```
conda activate verkko_v2.2.1
verkko -d run1_asm --hifi hifi.filt.fastq.gz --telomere-motif CCTAA
```
