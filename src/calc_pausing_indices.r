## File to calculate pausing indices using the method from Core, 2008

## Parse arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
    args <- c("--help")
}

## Help section
if("--help" %in% args) {
    cat("
      R Argument Parsing
 
      Arguments:
      --srr=SRR000000    - SRR Number
      --help              - print this text
 
      Example:
      ./script.r --srr=SRR123456 \n\n")
    
    q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))

if(is.null(argsL)) {
    print("You must provide an SRR")
    q()
}

.libPaths( c( .libPaths(), "/Users/zama8258/R/") )
.libPaths()
message("Established libpaths")

suppressMessages(library(tidyverse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(groHMM))

no_cores <- detectCores()
message(str_c("Using, ", no_cores, " cores."))

## Make sure we parse the correct srr
message(str_c("SRR is: ", argsL))


## Get refseq annotations if we don't have them
if(!file.exists("/scratch/Users/zama8258/hg19RefGene.sqlite")) {
    message("Downloading refseq annotations")
    rgdb <- makeTxDbFromUCSC(genome="hg19", tablename="ncbiRefSeq")
    saveDb(rgdb, file="/scratch/Users/zama8258/hg19RefGene.sqlite")
}

## Load refseq annotations if we do have them
message("Loading refseq annotations")
rgdb <- loadDb("/scratch/Users/zama8258/hg19RefGene.sqlite")

## Get transcripts from refseq
ts <- transcripts(rgdb)

## Parse our SRR bam file as a Granges object
message("Reading SRR BAM")
reads <- as(readGAlignments("/scratch/Shares/public/nascentdb/processedv2.0/bams/SRR2084576.trimmed.bam"), "GRanges")

## Calculate pausing indices using GROHmm code
message("Calculating Pausing Indices")
pi <- pausingIndex(ts, reads, size = 50, up = 1000, down = 1000, mc.cores=no_cores)

## Save for later analysis
## TODO - Export with SRR Name
message("Saving pausing index data")
save(pi, file="/scratch/Users/zama8258/pi_test.Rda")
