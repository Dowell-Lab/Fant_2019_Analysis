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
suppressMessages(library("argparse"))
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

df_sense <- countsSense %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_sense$coord = as.numeric(df_sense$coord)

df_antisense <- countsAntiSense %>% separate(Geneid, into = c("geneid", "coord"), sep = "/")
df_antisense$coord = as.numeric(df_antisense$coord)

## Normalize the counts for each region by TPM
tpm_normalize <- function(sense, antisense, sizefactor, sample) {
    ## First, calculate RPK for the strand only
    RPK_sample <- (sense[[sample]] / (sense$Length / (10 ^ 3)))
    ## Then, calculate RPK for both strands (for scaling factor)
    RPK <- ((sense[[sample]] + antisense[[sample]]) /
            (sense$Length / (10 ^ 3)))
    ## Then, calculate the scaling factor
    scale <- sum(RPK) / 1000000
    ## Divide RPK values by scaling factor
    out <- sizefactor * (RPK_sample / scale)
    return(out)
}

## TPM Normalize Sense Strand
df_sense$C413_1_S3_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, 1.0762, "C413_1_S3_R1_001.sorted.bam")
df_sense$C413_2_S4_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, 0.9514, "C413_2_S4_R1_001.sorted.bam")
df_sense$PO_1_S1_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, 0.7735, "PO_1_S1_R1_001.sorted.bam")
df_sense$PO_2_S2_R1_001.sorted.bam <-
    tpm_normalize(df_sense, df_antisense, 1.2678, "PO_2_S2_R1_001.sorted.bam")

## Repeat for Antisense
## df_antisense$C413_1_S3_R1_001.sorted.bam <-
##     tpm_normalize(df_antisense, df_sense, 1.0762, "C413_1_S3_R1_001.sorted.bam")
## df_antisense$C413_2_S4_R1_001.sorted.bam <-
##     tpm_normalize(df_antisense, df_sense, 0.9514, "C413_2_S4_R1_001.sorted.bam")
## df_antisense$PO_1_S1_R1_001.sorted.bam <-
##     tpm_normalize(df_antisense, df_sense, 0.7735, "PO_1_S1_R1_001.sorted.bam")
## df_antisense$PO_2_S2_R1_001.sorted.bam <-
##     tpm_normalize(df_antisense, df_sense, 1.2678, "PO_2_S2_R1_001.sorted.bam")

## Transform data by calculating the max, min, and mean for each row
df_sense <- df_sense %>%
    mutate(treat_mu = rowMeans(select(.,
                                      C413_1_S3_R1_001.sorted.bam,
                                      C413_2_S4_R1_001.sorted.bam)))
df_sense <- df_sense %>%
    mutate(control_mu = rowMeans(select(.,
                                        PO_1_S1_R1_001.sorted.bam,
                                        PO_2_S2_R1_001.sorted.bam)))
df_sense <- df_sense %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                            C413_1_S3_R1_001.sorted.bam,
                                            C413_2_S4_R1_001.sorted.bam,
                                            PO_1_S1_R1_001.sorted.bam,
                                PO_2_S2_R1_001.sorted.bam))
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

## library('fitdistrplus')
## treat_fit <- fitdist(final_sense$mean_treat_mu, 'gamma')
## treat_shape <- summary(treat_fit)$estimate[1]
## treat_rate <- summary(treat_fit)$estimate[2]
## control_fit <- fitdist(final_sense$mean_control_mu, 'gamma')
## control_shape <- summary(control_fit)$estimate[1]
## control_rate <- summary(control_fit)$estimate[2]
## final_sense <- final_sense %>%
##     mutate(mean_treat_mu = pgamma(mean_treat_mu, shape = treat_shape, rate = treat_rate))
## final_sense <- final_sense %>%
##     mutate(mean_control_mu = pgamma(mean_control_mu, shape = control_shape, rate = control_rate))

## df_antisense <- df_antisense %>%
##     mutate(treat_mu = rowMeans(select(.,
##                                       C413_1_S3_R1_001.sorted.bam,
##                                       C413_2_S4_R1_001.sorted.bam)))
## df_antisense <- df_antisense %>%
##     mutate(control_mu = rowMeans(select(.,
##                                         PO_1_S1_R1_001.sorted.bam,
##                                         PO_2_S2_R1_001.sorted.bam)))
## df_antisense <- df_antisense %>% subset(select = -c(geneid, Chr, Start, End, Length,
##                                                     C413_1_S3_R1_001.sorted.bam,
##                                             C413_2_S4_R1_001.sorted.bam,
##                                             PO_1_S1_R1_001.sorted.bam,
##                                             PO_2_S2_R1_001.sorted.bam))
## df_antisense <- df_antisense %>% mutate(coord = ifelse(Strand == '-', (num_bins - 1) - coord, coord))

## final_antisense <- df_antisense %>% group_by(coord) %>%
##     summarise(mean_treat_mu = -mean(treat_mu),
##               mean_treat_var = sd(treat_mu) / sqrt(length(treat_mu)),
##               mean_control_mu = -mean(control_mu),
##               mean_control_var = sd(control_mu) / sqrt(length(treat_mu))) %>%
##     mutate(mean_treat_max = mean_treat_mu + mean_treat_var,
##            mean_treat_min = mean_treat_mu - mean_treat_var,
##            mean_control_max = mean_control_mu + mean_control_var,
##            mean_control_min = mean_control_mu - mean_control_var)

## Set the correct axes based on filename
if (grepl('tss', sense)) {
    axes <- scale_x_continuous(breaks=c(0, 25, 50, 75, 100),
                               labels=c("-2kb", "-1kb", "TSS", "+1kb", "+2kb"),
                               limits=c(25,75))
    vert <- geom_vline(xintercept = 50, color = "black", linetype= "solid")
    ## pause <- geom_vline(xintercept = 51.5, color = "blue", linetype= "dashed", size = 0.4, alpha = 0.5)
    pause <- geom_blank()
} else if (grepl('body', sense)) {
    axes <- scale_x_continuous(breaks=c(0, 100),
                               labels=c("+2kb TSS", "-2kb TES"))
    vert <- geom_blank()
    pause <- geom_blank()
} else if (grepl('tes', sense)) {
    axes <- scale_x_continuous(breaks=c(0, 25, 50, 75, 100),
                               labels=c("-2kb", "-1kb", "TES", "+1kb", "+2kb"),
                               limits=c(25,75))
    vert <- geom_vline(xintercept = 50, color = "black", linetype= "dashed", size = 0.1, alpha=0.5)
    pause <- geom_blank()
} else {
    print("Pattern Matching Error In Plot Generation")
    axes <- geom_blank()
    vert <- geom_blank()
    pause <- geom_blank()
}

## Full Plot
library("ggthemes")
knockdownNorm <- factor(c('Knockdown', 'Control'))
ggplot() + theme_tufte() +
    scale_color_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    scale_fill_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    ## Sense Knockdown
    geom_line(data = final_sense, aes(x = coord,
                                      y = mean_treat_mu,
                                      color = 'Knockdown'), size = 0.5) +
    geom_ribbon(data = final_sense, aes(x = coord,
                                        ymin = mean_treat_min,
                                        ymax = mean_treat_max,
                                        fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
    ## Sense Control
    geom_line(data = final_sense, aes(x = coord,
                                      y = mean_control_mu,
                                      color = 'Control'), size = 0.5) +
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
    ## Plot line and axis adjustment
    axes + vert + pause + ylim(0, 4) +
    ## geom_hline(yintercept = 0, color = "black", size = 0.125) +
    ## Labels and such
    labs(title = "Metagene Plot", color = "Condition", fill = "Std. Dev. Mean") +
    xlab("Genomic Position") +
    ylab("Normalized Read Depth")

ggsave(outfile, width = 10, height = 5)

ggplot(data = final_sense) +
    geom_density(aes(mean_treat_mu)) +
    geom_density(aes(mean_control_mu)) +
    theme_tufte() +
    ggsave(str_c(outfile, "_dist.pdf"), width = 10, height = 5)

######################################################################
### metagene_graph_custom.r ends here
