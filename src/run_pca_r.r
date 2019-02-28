### run_pca_r.r --- run pca with batch correction
##
## Filename: run_pca_r.r
## Description: Run PCA
## Author: Student Zachary Maas <zama8258@colorado.edu>
## Maintainer: Student Zachary Maas <zama8258@colorado.edu>
## Created: Thu Jan 31 14:25:39 2019 (-0700)
##
######################################################################
##
### Commentary:
##
## Runs PCA with batch correction...
##
### Code:

library("tidyverse")
library("ggfortify")
library("ggplot2")
library("tools") # For debugging
library("limma")
library("sva")

setwd('/home/zach/dowell_lab/pausing_meta_analysis/out/counts')
counts <- read_delim("counts_top500_fullgene.txt", delim = "\t") %>%
    subset(select = -c(Geneid, Chr, Start, End, Strand, Length))

batch <- c("Day 1", "Day 2", "Day 1", "Day 2")
fx <- removeBatchEffect(counts, batch)
## fx <- ComBat(as.matrix(counts), batch)
counts_t <- data.frame(t(fx))
counts_t_filt <- counts_t[, which(apply(counts_t, 2, var) != 0)]

pca <- prcomp(counts_t_filt)

pcap <- data.frame(pca$x)
pcap$group <- c("Treatment", "Treatment", "Control", "Control")

ggplot(data = pcap, mapping = aes(x = PC1, y = PC2, color = group)) +
    geom_point() +
    labs(x = "PC1", y = "PC2", title = "Change in Pausing Index")
ggsave(str_c("/home/zach/dowell_lab/pausing_meta_analysis/out/prepublish/pca/pca_fullgene_top500.png"),
       plot = last_plot(),
       device = "png")
md5sum("/home/zach/dowell_lab/pausing_meta_analysis/out/prepublish/pca/pca_fullgene_top500.png")

######################################################################
### run_pca_r.r ends here
