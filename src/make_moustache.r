### make_moustache.r ---
##
## Filename: make_moustache.r
## Description:
## Author:
## Maintainer:
## Created: Fri Apr 26 15:31:13 2019 (-0600)
## Version:
## URL:
##
######################################################################
##
### Commentary:
##
## Make a moustache plot from gsea results
## 
### Code:

library('tidyverse')
library('ggthemes')
library('ggsci')

setwd('/home/zach/dowell_lab/pausing_meta_analysis/out/prepublish/gsea/output')

datp <- read_delim('gsea_report_for_na_pos_1548372740009.xls', delim = '\t')
datp$key <- ifelse(datp$`FDR q-val` < 0.1, "Significant Pos", "Not Significant")
datn <- read_delim('gsea_report_for_na_neg_1548372740009.xls', delim = '\t')
datn$key <- ifelse(datn$`FDR q-val` < 0.1, "Significant Neg", "Not Significant")
dat <- rbind(datp, datn)
dat$key <- as.factor(dat$key)

pal <- c("Significant Pos" = "red",
         "Significant Neg" = "green",
         "Not Significant" = "black")

# Anything under 0.1 significant, green on left, red on right
ggplot() +
    geom_point(data = dat, aes(x = `NES`, y = `FDR q-val`, color = key),
               alpha = 1, size = 4) +
    geom_hline(aes(yintercept = 0.1)) + geom_vline(aes(xintercept = 0)) +
    ## geom_point(data = dat, aes(x = `FDR q-val`, y = `NOM p-val`),
    ##            alpha = 0.5) +
    theme_tufte() +
    scale_color_manual(
        values = pal,
        limits = names(pal)) +
    labs(title = "Treatment vs Control",
         x = "Normalized Enrichment Score", y = "FDR q-value") +
    guides(color = FALSE) +
    ggsave('/home/zach/dowell_lab/pausing_meta_analysis/out/prepublish/gsea/moustache.png',
           device = 'png',
           width = 10,
           height = 5)

######################################################################
### make_moustache.r ends here
