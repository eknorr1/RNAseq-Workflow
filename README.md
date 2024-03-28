# RNAseq-Workflow

1. Download fastq files

  
2. Merge Reads 

3. Trim files with trim galore
trim_galore -q 20 --paired .fastq cont1_2.fastq
ls *_R1.fastq.gz | cut -d "_" -f 1 > seqlist
for filn in `cat seqlist`; do trim_galore -q 20 --paired $filn"_R1.fastq.gz" $filn"_R2.fastq.gz"; done

#bwa to map reads
bwa index pa14.fna 
for filn in `cat seqlist`; do bwa mem pa14.fna $filn"_R1_val_1.fq.gz" $filn"_R2_val_2.fq.gz" | samtools sort | samtools view -F 4 -o $filn".sorted.bam"; done
for filn in *sorted.bam; do samtools index $filn; done
grep -P "\tCDS\t" GCF_006974045.1_ASM697404v1_genomic.gff | cut -f 1,4,5 > pa14.bed
multiBamCov -bams 151.sorted.bam 221.sorted.bam 222.sorted.bam -bed pa14.bed > pa-bwa.counts.txt
grep -P "\tCDS\t" GCF_006974045.1_ASM697404v1_genomic.gff | cut -f 9 | cut -d "=" -f 3 | sed "s/gene-PA_//" | sed "s/;.*//" > pa_tags.txt
paste pa_tags.txt pa-bwa.counts.txt | cut -f 1,5-8 > pa-bwa.countsR.txt
sed -i 's/gene-EIP97_//' pa-bwa.countsR.txt
#need to nano `pa-bwa.countsR.txt` and add headers or export it to excel then do it and save as .csv

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
setwd("/Users/calebmallery/Desktop/")
countsTable = read.table("pa-bwa.countsR.csv", header=T, sep=",", row.names=1)
ColData = read.table("ColData.txt", header=T, sep="\t", row.names=1)
dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = ColData, design = ~ Condition)
dds <- dds[rowSums(counts(dds)) >100, ]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Condition", "ko", "control"))
write.table(res, file="151_DESeq2.txt", sep='\t')
