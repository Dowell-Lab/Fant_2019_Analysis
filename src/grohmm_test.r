## source("http://bioconductor.org/biocLite.R")
## biocLite("groHMM")
.libPaths( c( .libPaths(), "/Users/zama8258/R/") )
.libPaths()
message("Established libpaths")

library(tidyverse)

library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(groHMM)

no_cores <- detectCores()

## For refseq annotations:
if(!file.exists("/scratch/Users/zama8258/hg19RefGene.sqlite")) {
    rgdb <- makeTxDbFromUCSC(genome="hg19", tablename="ncbiRefSeq")
    saveDb(rgdb, file="/scratch/Users/zama8258/hg19RefGene.sqlite")
}

rgdb <- loadDb("/scratch/Users/zama8258/hg19RefGene.sqlite")

ts <- transcripts(rgdb)

reads <- as(readGAlignments("/scratch/Shares/public/nascentdb/processedv2.0/bams/SRR2084576.trimmed.bam"), "GRanges")

pi <- pausingIndex(ts, reads, size = 50, up = 1000, down = 1000, mc.cores=no_cores)

save(pi, file="/scratch/Users/zama8258/pi_test.Rda")
