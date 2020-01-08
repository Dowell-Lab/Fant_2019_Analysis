### metagene_graph_custom.r --- Custom Metagene Plots
##
## Filename: metagene_graph_custom.r
## Author: Student Zachary Maas <zama8258@colorado.edu>
## Maintainer: Student Zachary Maas <zama8258@colorado.edu>
##
######################################################################
##
### Commentary:
##
## Builds Metagene Plots from Custom CountsSense File
##
######################################################################
##
### Code:

suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
## suppressMessages(library("argparse"))
suppressMessages(library("matrixStats"))

parser <- ArgumentParser()
parser$add_argument("-s", "--sense", action="store", dest="countsSense",
                    help="The sense reads counts table to use")
parser$add_argument("-a", "--antisense", action="store", dest="countsAntiSense",
                    help="The antisense reads table to use")
parser$add_argument("-n", "--numbins", action="store", dest="numbins",
                    help="The number of bins used")
parser$add_argument("-o", "--output", action="store", dest="outfile",
                    help="The output image file name.")
args <- parser$parse_args()

sense <- args$countsSense
antisense <- args$countsAntiSense
num_bins <- as.numeric(args$numbins)
outfile <- args$outfile

countsSense <- read_delim(sense, delim="\t")
countsAntiSense <- read_delim(antisense, delim="\t")

setwd('/home/zach/dowell_lab/pausing_meta_analysis/out/counts/drosophila')
countsSense <- read_delim('metagene_counts_sense_fix.txt', delim="\t")
countsAntiSense <- read_delim('metagene_counts_antisense_fix.txt', delim="\t")
outfile <- '/home/zach/dowell_lab/pausing_meta_analysis/out/drosophila/testing.pdf'
num_bins <- 25

df_sense <- countsSense %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_sense$coord = as.numeric(df_sense$coord)

df_antisense <- countsAntiSense %>% separate(Geneid, into = c("geneid", "coord"), sep = "/")
df_antisense$coord = as.numeric(df_antisense$coord)

## Normalize the counts for each region by TPM
tpm_normalize <- function(sense, antisense, sample) {
    ## First, calculate RPK for the strand only
    RPK_sample <- (sense[[sample]] / (sense$Length / (10 ^ 3)))
    ## Then, calculate RPK for both strands (for scaling factor)
    RPK <- ((sense[[sample]] + antisense[[sample]]) /
            (sense$Length / (10 ^ 3)))
    ## Then, calculate the scaling factor
    scale <- sum(RPK) / 1000000
    ## Divide RPK values by scaling factor
    out <- RPK_sample / scale
    return(out)
}

## TPM Normalize Sense Strand
df_sense$Control_1_S1_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, "Control_1_S1_R1_001.sorted.bam")
df_sense$Control_2_S2_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, "Control_2_S2_R1_001.sorted.bam")
df_sense$Control_3_S3_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, "Control_3_S3_R1_001.sorted.bam")
df_sense$Taf_1_S4_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, "Taf_1_S4_R1_001.sorted.bam")
df_sense$Taf_2_S5_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, "Taf_2_S5_R1_001.sorted.bam")
df_sense$Taf_3_S6_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, "Taf_3_S6_R1_001.sorted.bam")

##  FIX TYPO
## upregulated <- df_sense %>%
##     group_by(geneid) %>%
##     summarise(c1 = sum(Control_1_S1_R1_001.sorted.bam),
##               c2 = sum(Control_2_S2_R1_001.sorted.bam),
##               c3 = sum(Control_3_S3_R1_001.sorted.bam),
##               p1 = sum(Taf_1_S4_R1_001.sorted.bam),
##               p2 = sum(Taf_2_S5_R1_001.sorted.bam),
##               p3 = sum(Taf_3_S6_R1_001.sorted.bam)) %>%
##     mutate(cmean = rowMeans(select(., c1, c2, c3)),
##            pmean = rowMeans(select(., p1, p2, p3))) %>%
##     mutate(expressionchange = pmean - cmean) %>%
##     subset(select = -c(c1, c2, c3, p1, p2, p3, cmean, pmean)) %>%
##     filter(expressionchange > 0)
## df_sense <- left_join(upregulated, df_sense)

## Repeat for Antisense
df_antisense$Control_1_S1_R1_001.sorted.bam <-
    tpm_normalize(df_antisense, df_sense, "Control_1_S1_R1_001.sorted.bam")
df_antisense$Control_2_S2_R1_001.sorted.bam <-
    tpm_normalize(df_antisense, df_sense, "Control_2_S2_R1_001.sorted.bam")
df_antisense$Control_3_S3_R1_001.sorted.bam <-
    tpm_normalize(df_antisense, df_sense, "Control_3_S3_R1_001.sorted.bam")
df_antisense$Taf_1_S4_R1_001.sorted.bam <-
    tpm_normalize(df_antisense, df_sense, "Taf_1_S4_R1_001.sorted.bam")
df_antisense$Taf_2_S5_R1_001.sorted.bam <-
    tpm_normalize(df_antisense, df_sense, "Taf_2_S5_R1_001.sorted.bam")
df_antisense$Taf_3_S6_R1_001.sorted.bam <-
    tpm_normalize(df_antisense, df_sense, "Taf_3_S6_R1_001.sorted.bam")

## Transform data by calculating the max, min, and mean for each row
df_sense <- df_sense %>%
    mutate(control_mu = rowMeans(select(.,
                                        Control_1_S1_R1_001.sorted.bam,
                                      Control_2_S2_R1_001.sorted.bam,
                                      Control_3_S3_R1_001.sorted.bam)))
df_sense <- df_sense %>%
    mutate(treat_mu = rowMeans(select(.,
                                      Taf_1_S4_R1_001.sorted.bam,
                                        Taf_2_S5_R1_001.sorted.bam,
                                        Taf_3_S6_R1_001.sorted.bam)))
df_sense <- df_sense %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                            Control_1_S1_R1_001.sorted.bam,
                                            Control_2_S2_R1_001.sorted.bam,
                                            Control_3_S3_R1_001.sorted.bam,
                                            Taf_1_S4_R1_001.sorted.bam,
                                            Taf_2_S5_R1_001.sorted.bam,
                                            Taf_3_S6_R1_001.sorted.bam))
df_sense <- df_sense %>% mutate(coord = ifelse(Strand == '-', (num_bins - 1) - coord, coord))

final_sense <- df_sense %>% group_by(coord) %>%
    summarise(mean_treat_mu = mean(treat_mu),
              mean_treat_var = sd(treat_mu) / sqrt(length(treat_mu)),
              mean_control_mu = mean(control_mu),
              mean_control_var = sd(control_mu) / sqrt(length(treat_mu))) %>%
    mutate(mean_treat_max = mean_treat_mu + mean_treat_var,
           mean_treat_min = mean_treat_mu - mean_treat_var,
           mean_control_max = mean_control_mu + mean_control_var,
           mean_control_min = mean_control_mu - mean_control_var)


df_antisense <- df_antisense %>%
    mutate(control_mu = rowMeans(select(.,
                                        Control_1_S1_R1_001.sorted.bam,
                                      Control_2_S2_R1_001.sorted.bam,
                                      Control_3_S3_R1_001.sorted.bam)))
df_antisense <- df_antisense %>%
    mutate(treat_mu = rowMeans(select(.,
                                      Taf_1_S4_R1_001.sorted.bam,
                                        Taf_2_S5_R1_001.sorted.bam,
                                        Taf_3_S6_R1_001.sorted.bam)))
df_antisense <- df_antisense %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                                    Control_1_S1_R1_001.sorted.bam,
                                                    Control_2_S2_R1_001.sorted.bam,
                                                    Control_3_S3_R1_001.sorted.bam,
                                                    Taf_1_S4_R1_001.sorted.bam,
                                                    Taf_2_S5_R1_001.sorted.bam,
                                                    Taf_3_S6_R1_001.sorted.bam))
df_antisense <- df_antisense %>% mutate(coord = ifelse(Strand == '-', (num_bins - 1) - coord, coord))

final_antisense <- df_antisense %>% group_by(coord) %>%
    summarise(mean_treat_mu = -mean(treat_mu),
              mean_treat_var = sd(treat_mu) / sqrt(length(treat_mu)),
              mean_control_mu = -mean(control_mu),
              mean_control_var = sd(control_mu) / sqrt(length(treat_mu))) %>%
    mutate(mean_treat_max = mean_treat_mu + mean_treat_var,
           mean_treat_min = mean_treat_mu - mean_treat_var,
           mean_control_max = mean_control_mu + mean_control_var,
           mean_control_min = mean_control_mu - mean_control_var)

## Full Plot
library("ggthemes")
knockdownNorm <- factor(c('Knockdown', 'Control'))
ggplot() + theme_tufte() +
    scale_color_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    scale_fill_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    ## Sense Knockdown
    geom_line(data = final_sense, aes(x = coord,
                                      y = mean_treat_mu,
                                      color = 'Knockdown'), size = 0.25) +
    geom_ribbon(data = final_sense, aes(x = coord,
                                        ymin = mean_treat_min,
                                        ymax = mean_treat_max,
                                        fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
    ## Sense Control
    geom_line(data = final_sense, aes(x = coord,
                                      y = mean_control_mu,
                                      color = 'Control'), size = 0.25) +
    geom_ribbon(data = final_sense, aes(x = coord,
                                        ymin = mean_control_min,
                                        ymax = mean_control_max,
                                        fill = 'Control'), size=0.25, alpha = 0.2) +
    ## Antisense Knockdown
    ## geom_line(data = final_antisense, aes(x = coord,
    ##                                       y = mean_treat_mu,
    ##                                       color = 'Knockdown'), size=0.25) +
    ## geom_ribbon(data = final_antisense, aes(x = coord,
    ##                                         ymin = mean_treat_min,
    ##                                         ymax = mean_treat_max,
    ##                                         fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
    ## ## Antisense Control
    ## geom_line(data = final_antisense, aes(x = coord,
    ##                                       y = mean_control_mu,
    ##                                       color = 'Control'), size=0.25) +
    ## geom_ribbon(data = final_antisense, aes(x = coord,
    ##                                         ymin = mean_control_min,
    ##                                         ymax = mean_control_max,
    ##                                         fill = 'Control'), size=0.25, alpha = 0.2) +
    geom_hline(yintercept = 0, color = "black", size = 0.125) +
    labs(title = "Metagene Plot", color = "Condition", fill = "Std. Dev. Mean") +
    xlab("Bins") + ylab("Normalized Read Depth")

ggsave(outfile, width = 10, height = 5)

## ## Filter data for inset
## inset_sense <- final_sense %>% filter(coord > 40) %>% filter(coord < 80)
## inset_antisense <- final_antisense %>% filter(coord > 40) %>% filter(coord < 80)

## ## Inset Plot
## ggplot() + theme_tufte() +
##     scale_color_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
##     scale_fill_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
##     ## Sense Knockdown
##     geom_line(data = inset_sense, aes(x = coord,
##                                       y = mean_treat_mu,
##                                       color = 'Knockdown'), size = 0.25) +
##     geom_ribbon(data = inset_sense, aes(x = coord,
##                                         ymin = mean_treat_min,
##                                         ymax = mean_treat_max,
##                                         fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
##     ## Sense Control
##     geom_line(data = inset_sense, aes(x = coord,
##                                       y = mean_control_mu,
##                                       color = 'Control'), size = 0.25) +
##     geom_ribbon(data = inset_sense, aes(x = coord,
##                                         ymin = mean_control_min,
##                                         ymax = mean_control_max,
##                                         fill = 'Control'), size=0.25, alpha = 0.2) +
##     ## Antisense Knockdown
##     geom_line(data = inset_antisense, aes(x = coord,
##                                           y = mean_treat_mu,
##                                           color = 'Knockdown'), size=0.25) +
##     geom_ribbon(data = inset_antisense, aes(x = coord,
##                                             ymin = mean_treat_min,
##                                             ymax = mean_treat_max,
##                                             fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
##     ## Antisense Control
##     geom_line(data = inset_antisense, aes(x = coord,
##                                           y = mean_control_mu,
##                                           color = 'Control'), size=0.25) +
##     geom_ribbon(data = inset_antisense, aes(x = coord,
##                                             ymin = mean_control_min,
##                                             ymax = mean_control_max,
##                                             fill = 'Control'), size=0.25, alpha = 0.2) +
##     geom_hline(yintercept = 0, color = "black", size = 0.125) +
##     labs(title = "Metagene Plot", color = "Condition", fill = "Std. Dev. Mean") +
##     xlab("Bins") + ylab("Normalized Read Depth")

## ggsave(paste0(outfile, "_inset.pdf"), width = 10, height = 5)

######################################################################
### metagene_graph_custom.r ends here
