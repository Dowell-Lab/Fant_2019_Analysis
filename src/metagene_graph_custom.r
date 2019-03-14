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

## countsSense <- read_delim('metagene_counts_sense_fix.txt', delim="\t")
## countsAntiSense <- read_delim('metagene_counts_antisense_fix.txt', delim="\t")

df_sense <- countsSense %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_sense$coord = as.numeric(df_sense$coord)
df_sense <- df_sense %>%
    mutate(treat_mu = rowMeans(select(.,
                                      C413_1_S3_R1_001.sorted.bam,
                                      C413_2_S4_R1_001.sorted.bam)))
df_sense <- df_sense %>%
    mutate(treat_max = pmax(C413_1_S3_R1_001.sorted.bam,
                            C413_2_S4_R1_001.sorted.bam))
df_sense <- df_sense %>%
    mutate(treat_min = pmin(C413_1_S3_R1_001.sorted.bam,
                            C413_2_S4_R1_001.sorted.bam))
df_sense <- df_sense %>%
    mutate(control_mu = rowMeans(select(.,
                                        PO_1_S1_R1_001.sorted.bam,
                                        PO_2_S2_R1_001.sorted.bam)))
df_sense <- df_sense %>%
    mutate(control_max = pmax(PO_1_S1_R1_001.sorted.bam,
                              PO_2_S2_R1_001.sorted.bam))
df_sense <- df_sense %>%
    mutate(control_min = pmin(PO_1_S1_R1_001.sorted.bam,
                              PO_2_S2_R1_001.sorted.bam))
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

df_antisense <- countsAntiSense %>% separate(Geneid, into = c("geneid", "coord"), sep = "/")
df_antisense$coord = as.numeric(df_antisense$coord)
df_antisense <- df_antisense %>%
    mutate(treat_mu = rowMeans(select(.,
                                      C413_1_S3_R1_001.sorted.bam,
                                      C413_2_S4_R1_001.sorted.bam)))
df_antisense <- df_antisense %>%
    mutate(treat_max = pmax(C413_1_S3_R1_001.sorted.bam,
                            C413_2_S4_R1_001.sorted.bam))
df_antisense <- df_antisense %>%
    mutate(treat_min = pmin(C413_1_S3_R1_001.sorted.bam,
                            C413_2_S4_R1_001.sorted.bam))
df_antisense <- df_antisense %>%
    mutate(control_mu = rowMeans(select(.,
                                        PO_1_S1_R1_001.sorted.bam,
                                        PO_2_S2_R1_001.sorted.bam)))
df_antisense <- df_antisense %>%
    mutate(control_max = pmax(PO_1_S1_R1_001.sorted.bam,
                              PO_2_S2_R1_001.sorted.bam))
df_antisense <- df_antisense %>%
    mutate(control_min = pmin(PO_1_S1_R1_001.sorted.bam,
                              PO_2_S2_R1_001.sorted.bam))
df_antisense <- df_antisense %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                                    C413_1_S3_R1_001.sorted.bam,
                                            C413_2_S4_R1_001.sorted.bam,
                                            PO_1_S1_R1_001.sorted.bam,
                                            PO_2_S2_R1_001.sorted.bam))
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
ggplot() + theme_tufte() +
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
    geom_line(data = final_antisense, aes(x = coord,
                                          y = mean_treat_mu,
                                          color = 'Knockdown'), size=0.25) +
    geom_ribbon(data = final_antisense, aes(x = coord,
                                            ymin = mean_treat_min,
                                            ymax = mean_treat_max,
                                            fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
    ## Antisense Control
    geom_line(data = final_antisense, aes(x = coord,
                                          y = mean_control_mu,
                                          color = 'Control'), size=0.25) +
    geom_ribbon(data = final_antisense, aes(x = coord,
                                            ymin = mean_control_min,
                                            ymax = mean_control_max,
                                            fill = 'Control'), size=0.25, alpha = 0.2) +
    geom_hline(yintercept = 0, color = "black", size = 0.125) +
    labs(title = "Metagene Plot", color = "Condition", fill = "Std. Dev. Mean") +
    xlab("Bins") + ylab("Normalized Read Depth")
## ggsave("img.png", width = 10, height = 5)

ggsave(outfile, width = 10, height = 5)

## Filter data for inset
inset_sense <- final_sense %>% filter(coord > 40) %>% filter(coord < 80)
inset_antisense <- final_antisense %>% filter(coord > 40) %>% filter(coord < 80)

## Inset Plot
ggplot() + theme_tufte() +
    ## Sense Knockdown
    geom_line(data = inset_sense, aes(x = coord,
                                      y = mean_treat_mu,
                                      color = 'Knockdown'), size = 0.25) +
    geom_ribbon(data = inset_sense, aes(x = coord,
                                        ymin = mean_treat_min,
                                        ymax = mean_treat_max,
                                        fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
    ## Sense Control
    geom_line(data = inset_sense, aes(x = coord,
                                      y = mean_control_mu,
                                      color = 'Control'), size = 0.25) +
    geom_ribbon(data = inset_sense, aes(x = coord,
                                        ymin = mean_control_min,
                                        ymax = mean_control_max,
                                        fill = 'Control'), size=0.25, alpha = 0.2) +
    ## Antisense Knockdown
    geom_line(data = inset_antisense, aes(x = coord,
                                          y = mean_treat_mu,
                                          color = 'Knockdown'), size=0.25) +
    geom_ribbon(data = inset_antisense, aes(x = coord,
                                            ymin = mean_treat_min,
                                            ymax = mean_treat_max,
                                            fill = 'Knockdown'), size = 0.25, alpha = 0.2) +
    ## Antisense Control
    geom_line(data = inset_antisense, aes(x = coord,
                                          y = mean_control_mu,
                                          color = 'Control'), size=0.25) +
    geom_ribbon(data = inset_antisense, aes(x = coord,
                                            ymin = mean_control_min,
                                            ymax = mean_control_max,
                                            fill = 'Control'), size=0.25, alpha = 0.2) +
    geom_hline(yintercept = 0, color = "black", size = 0.125) +
    labs(title = "Metagene Plot", color = "Condition", fill = "Std. Dev. Mean") +
    xlab("Bins") + ylab("Normalized Read Depth")

ggsave(paste0("inset_", outfile), width = 10, height = 5)

######################################################################
### metagene_graph_custom.r ends here
