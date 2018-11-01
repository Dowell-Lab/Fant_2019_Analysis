## load library for genomic annotations
library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)

## load the transcript annotation file from UCSC.  Make sure to enter the correct genome version
txdb=makeTxDbFromUCSC(genome='hg38',tablename='refGene')

## Use the function transcriptsBy(txdb,'gene') for the whole genic region instead of just exons
ex_by_gene=exonsBy(txdb,'gene')

## load the samtools library for R
library(Rsubread)

## read the sequencing read alignment into R (combine with next step to save memory)
features <- featureCounts(c("C413_1_S3_R1_001.trim.sorted.bam",
                            "C413_2_S4_R1_001.trim.sorted.bam",
    "PO_1_S1_R1_001.trim.sorted.bam",
    "PO_2_S2_R1_001.trim.sorted.bam"), nthreads=16)

df <- as.tibble(features$counts)

## set the gene IDs to the table row names
rownames(countTable)=names(ex_by_gene)

##output tag counts to a file
write.table(countTable,file="countTable.txt",sep="\t")

##removing rows that are zero for all genes (edgeR and DESeq have trouble with these)
## x <- rowSums(countTable==0)!=ncol(countTable)
## newCountTable <- countTable[x,]
