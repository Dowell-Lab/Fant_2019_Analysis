message("Generating Figures")
suppressMessages(library("tidyverse"))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))

## Parse arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if (length(args) < 1) {
    args <- c("--help")
}

## Help section
if ("--help" %in% args) {
    cat("
Generate Figures for a Given SRR

Arguments:
   --srr=SRR000000 - SRR Number
   --help          - Print this help message
 
Example:
   ./calc_pausing_indices.r --srr=SRR123456 \n\n")
    
    q(save = "no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
srr <- argsL

if (is.null(srr)) {
    warningr("You must provide an SRR")
    q()
}

srr <- "SRR2084576"
message(str_c("SRR is: ", srr))

## Generate figures...

load(str_c("/scratch/Users/zama8258/pause_output/Pause_Core_", srr, ".Rda"))

## TODO Document these transformations better
message("Transforming Data")
dd <- as.tibble(pi) %>% filter(GeneID != "")

dd$GeneID = as.integer(as.character(dd$GeneID))

## Load refseq annotations if we do have them
message("Loading refseq annotations")
rgdb <- loadDb("/scratch/Users/zama8258/hg19RefGene.sqlite")

## Get transcripts from refseq
ts <- transcripts(rgdb)
tsd <- as.tibble(ts)
colnames(tsd)[colnames(tsd) == "tx_id"] <- "GeneID"

prd2000 <- read_delim(str_c("/scratch/Users/zama8258/pause_output/", srr, "_pause_ratios_2000.data"), 
    col_names = c("tx_name", "strand", "Pause_2000"), delim = " ")
prd5000 <- read_delim(str_c("/scratch/Users/zama8258/pause_output/", srr, "_pause_ratios_5000.data"), 
    col_names = c("tx_name", "strand", "Pause_5000"), delim = " ")

ddt <- as.tibble(merge(x = dd, y = tsd[, c("tx_name", "GeneID")], by.x = "GeneID", 
    by.y = "GeneID", all.x = TRUE))
ddt <- ddt[, c("GeneID", "tx_name", "Fisher", "Pause")]
ddt$tx_name <- substr(ddt$tx_name, 1, nchar(ddt$tx_name) - 2)
ddt <- ddt %>% arrange(tx_name)
prd2000 <- prd2000 %>% arrange(tx_name)
prd5000 <- prd5000 %>% arrange(tx_name)
ddt <- left_join(x = ddt, y = prd2000[, c("tx_name", "Pause_2000")], by = "tx_name")
ddt <- left_join(x = ddt, y = prd5000[, c("tx_name", "Pause_5000")], by = "tx_name")
ddt <- na.omit(ddt)

## ddt <- ddt %>% filter(Fisher < 0.05)

ddf <- ddt %>% filter(Pause < 200, Pause_5000 < 200)

library(ggthemes)

## Figures to Generate

## 2k vs 5k
ggplot(data = ddt, mapping = aes(x = Pause_2000, y = Pause_5000), alpha = 1/10) + 
    geom_jitter() + geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + 
    labs(x = "2kb Fixed Window", y = "5kb Fixed Window", title = "Comparison of Pausing Index Methods")
ggsave(str_c("/scratch/Users/zama8258/pause_output/", srr, "_comparison_2k_5k.png"), 
    plot = last_plot(), device = "png")

## Core vs 5k
ggplot(data = ddf, mapping = aes(x = Pause, y = Pause_5000), alpha = 1/10) + geom_jitter() + 
    geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + labs(x = "Core 2008 Method", 
    y = "5kb Fixed Window", title = "Comparison of Pausing Index Methods")
ggsave(str_c("/scratch/Users/zama8258/pause_output/", srr, "_comparison_Core_5k.png"), 
    plot = last_plot(), device = "png")

## Core vs 2k
ggplot(data = ddt, mapping = aes(x = Pause, y = Pause_2000), alpha = 1/10) + geom_jitter() + 
    geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + labs(x = "Core 2008 Method", 
    y = "2kb Fixed Window", title = "Comparison of Pausing Index Methods")
ggsave(str_c("/scratch/Users/zama8258/pause_output/", srr, "_comparison_Core_2k.png"), 
    plot = last_plot(), device = "png")
