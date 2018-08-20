## File to calculate pausing indices using TFit annotations

## Parse arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
    args <- c("--help")
}

## Help section
if("--help" %in% args) {
    cat("
Calculate Pausing Indices using Core 2008 Method

Arguments:
   --srr=SRR000000 - SRR Number
   --help          - Print this help message

Example:
   ./calc_pausing_indices.r --srr=SRR123456 \n\n")

    q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
srr <- argsL

if(is.null(srr)) {
    warningr("You must provide an SRR")
    q()
}

message("Libpaths:")
.libPaths( c( .libPaths(), "/Users/zama8258/R/") )
.libPaths()
message("Established libpaths")

suppressMessages(library(tidyverse))
## suppressMessages(library(GenomicFeatures))
## suppressMessages(library(GenomicRanges))
## suppressMessages(library(rtracklayer))
## suppressMessages(library(groHMM))
library(tfitParser)

no_cores <- detectCores()
message(str_c("Using, ", no_cores, " cores."))

## Make sure we parse the correct srr
message(str_c("SRR is: ", srr))

tf <- tfitParser::parse_outputfile("/scratch/Users/zama8258/tfit_runs/tfit_test/model/SRR2084576.tfit_mod_old9-1_K_models_MLE.tsv")

single <- tf[[2]]
single$Mu = as.integer(single$Mu)
single$Sigma = as.double(single$Sigma)
single$Lambda = as.double(single$Lambda)
single$Pi = as.double(single$Pi)
single$Fp = as.integer(single$Fp)
single$W_fk = as.double(single$W_fk)
single$W_rk = as.double(single$W_rk)
single$B = as.integer(single$B)
single$A = as.integer(single$A)

single <- single %>% filter(Log_Likelihood >= 0)
