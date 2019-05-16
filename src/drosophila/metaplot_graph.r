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
## Builds Metagene Plots from Custom CountsProximal File
##
######################################################################
##
### Code:

suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("matrixStats"))

## parser <- ArgumentParser()
## parser$add_argument("-s", "--proximal", action="store", dest="countsProximal",
##                     help="The proximal reads counts table to use")
## parser$add_argument("-a", "--distal", action="store", dest="countsDistal",
##                     help="The distal reads table to use")
## parser$add_argument("-n", "--numbins", action="store", dest="numbins",
##                     help="The number of bins used")
## parser$add_argument("-o", "--output", action="store", dest="outfile",
##                     help="The output image file name.")
## args <- parser$parse_args()

## proximal <- args$countsProximal
## distal <- args$countsDistal
## num_bins <- as.numeric(args$numbins)
## outfile <- args$outfile

## countsProximal <- read_delim(proximal, delim="\t")
## countsDistal <- read_delim(distal, delim="\t")

setwd("/home/zach/dowell_lab/pausing_meta_analysis/out/drosophila/metaplot")
num_bins=140
countsTop <- read_delim('topGenes_metagene_counts_sense_fix.txt', delim="\t")
countsPaused <- read_delim('genelist-paused_metagene_counts_sense_fix.txt', delim="\t")
countsProximal <- read_delim('genelist-Prox_metagene_counts_sense_fix.txt', delim="\t")
countsDistal <- read_delim('genelist-Dist_metagene_counts_sense_fix.txt', delim="\t")

df_top <- countsTop %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_top$coord = as.numeric(df_top$coord)

df_paused <- countsPaused %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_paused$coord = as.numeric(df_paused$coord)

df_proximal <- countsProximal %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_proximal$coord = as.numeric(df_proximal$coord)

df_distal <- countsDistal %>% separate(Geneid, into = c("geneid", "coord"), sep = "/")
df_distal$coord = as.numeric(df_distal$coord)

## Normalize the counts for each region by TPM
tpm_normalize <- function(proximal, sample) {
    ## First, calculate RPK for the strand only
    RPK_sample <- (proximal[[sample]] / (proximal$Length / (10 ^ 3)))
    ## Then, calculate RPK for both strands (for scaling factor)
    RPK <- ((proximal[[sample]]) /
            (proximal$Length / (10 ^ 3)))
    ## Then, calculate the scaling factor
    scale <- sum(RPK) / 1000000
    ## Divide RPK values by scaling factor
    out <- RPK_sample / scale
    return(out)
}

df_top$Control_1_S1_R1_001.sorted.bam <-
    tpm_normalize(df_top, "Control_1_S1_R1_001.sorted.bam")
df_top$Control_2_S2_R1_001.sorted.bam <-
    tpm_normalize(df_top, "Control_2_S2_R1_001.sorted.bam")
df_top$Control_3_S3_R1_001.sorted.bam <-
    tpm_normalize(df_top, "Control_3_S3_R1_001.sorted.bam")
df_top$Taf_1_S4_R1_001.sorted.bam <-
    tpm_normalize(df_top, "Taf_1_S4_R1_001.sorted.bam")
df_top$Taf_2_S5_R1_001.sorted.bam <-
    tpm_normalize(df_top, "Taf_2_S5_R1_001.sorted.bam")
df_top$Taf_3_S6_R1_001.sorted.bam <-
    tpm_normalize(df_top, "Taf_3_S6_R1_001.sorted.bam")

df_paused$Control_1_S1_R1_001.sorted.bam <-
    tpm_normalize(df_paused, "Control_1_S1_R1_001.sorted.bam")
df_paused$Control_2_S2_R1_001.sorted.bam <-
    tpm_normalize(df_paused, "Control_2_S2_R1_001.sorted.bam")
df_paused$Control_3_S3_R1_001.sorted.bam <-
    tpm_normalize(df_paused, "Control_3_S3_R1_001.sorted.bam")
df_paused$Taf_1_S4_R1_001.sorted.bam <-
    tpm_normalize(df_paused, "Taf_1_S4_R1_001.sorted.bam")
df_paused$Taf_2_S5_R1_001.sorted.bam <-
    tpm_normalize(df_paused, "Taf_2_S5_R1_001.sorted.bam")
df_paused$Taf_3_S6_R1_001.sorted.bam <-
    tpm_normalize(df_paused, "Taf_3_S6_R1_001.sorted.bam")

## TPM Normalize Proximal
df_proximal$Control_1_S1_R1_001.sorted.bam <-
    tpm_normalize(df_proximal, "Control_1_S1_R1_001.sorted.bam")
df_proximal$Control_2_S2_R1_001.sorted.bam <-
    tpm_normalize(df_proximal, "Control_2_S2_R1_001.sorted.bam")
df_proximal$Control_3_S3_R1_001.sorted.bam <-
    tpm_normalize(df_proximal, "Control_3_S3_R1_001.sorted.bam")
df_proximal$Taf_1_S4_R1_001.sorted.bam <-
    tpm_normalize(df_proximal, "Taf_1_S4_R1_001.sorted.bam")
df_proximal$Taf_2_S5_R1_001.sorted.bam <-
    tpm_normalize(df_proximal, "Taf_2_S5_R1_001.sorted.bam")
df_proximal$Taf_3_S6_R1_001.sorted.bam <-
    tpm_normalize(df_proximal, "Taf_3_S6_R1_001.sorted.bam")

## Repeat for Distal
df_distal$Control_1_S1_R1_001.sorted.bam <-
    tpm_normalize(df_distal, "Control_1_S1_R1_001.sorted.bam")
df_distal$Control_2_S2_R1_001.sorted.bam <-
    tpm_normalize(df_distal, "Control_2_S2_R1_001.sorted.bam")
df_distal$Control_3_S3_R1_001.sorted.bam <-
    tpm_normalize(df_distal, "Control_3_S3_R1_001.sorted.bam")
df_distal$Taf_1_S4_R1_001.sorted.bam <-
    tpm_normalize(df_distal, "Taf_1_S4_R1_001.sorted.bam")
df_distal$Taf_2_S5_R1_001.sorted.bam <-
    tpm_normalize(df_distal, "Taf_2_S5_R1_001.sorted.bam")
df_distal$Taf_3_S6_R1_001.sorted.bam <-
    tpm_normalize(df_distal, "Taf_3_S6_R1_001.sorted.bam")

## Transform data by calculating the max, min, and mean for each row
df_top <- df_top %>%
    mutate(treat_mu = rowMeans(select(.,
                                      Control_1_S1_R1_001.sorted.bam,
                                      Control_2_S2_R1_001.sorted.bam,
                                      Control_3_S3_R1_001.sorted.bam)))
df_top <- df_top %>%
    mutate(control_mu = rowMeans(select(.,
                                        Taf_1_S4_R1_001.sorted.bam,
                                        Taf_2_S5_R1_001.sorted.bam,
                                        Taf_3_S6_R1_001.sorted.bam)))
df_top <- df_top %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                        Control_1_S1_R1_001.sorted.bam,
                                              Control_2_S2_R1_001.sorted.bam,
                                              Control_3_S3_R1_001.sorted.bam,
                                              Taf_1_S4_R1_001.sorted.bam,
                                              Taf_2_S5_R1_001.sorted.bam,
                                              Taf_3_S6_R1_001.sorted.bam))
final_top <- df_top %>% group_by(coord) %>%
    summarise(mean_treat_mu = mean(treat_mu),
              mean_treat_var = sd(treat_mu) / sqrt(length(treat_mu)),
              mean_control_mu = mean(control_mu),
              mean_control_var = sd(control_mu) / sqrt(length(treat_mu))) %>%
    mutate(mean_treat_max = mean_treat_mu + mean_treat_var,
           mean_treat_min = mean_treat_mu - mean_treat_var,
           mean_control_max = mean_control_mu + mean_control_var,
           mean_control_min = mean_control_mu - mean_control_var)

df_paused <- df_paused %>%
    mutate(treat_mu = rowMeans(select(.,
                                      Control_1_S1_R1_001.sorted.bam,
                                      Control_2_S2_R1_001.sorted.bam,
                                      Control_3_S3_R1_001.sorted.bam)))
df_paused <- df_paused %>%
    mutate(control_mu = rowMeans(select(.,
                                        Taf_1_S4_R1_001.sorted.bam,
                                        Taf_2_S5_R1_001.sorted.bam,
                                        Taf_3_S6_R1_001.sorted.bam)))
df_paused <- df_paused %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                              Control_1_S1_R1_001.sorted.bam,
                                                  Control_2_S2_R1_001.sorted.bam,
                                                  Control_3_S3_R1_001.sorted.bam,
                                                  Taf_1_S4_R1_001.sorted.bam,
                                                  Taf_2_S5_R1_001.sorted.bam,
                                                  Taf_3_S6_R1_001.sorted.bam))
final_paused <- df_paused %>% group_by(coord) %>%
    summarise(mean_treat_mu = mean(treat_mu),
              mean_treat_var = sd(treat_mu) / sqrt(length(treat_mu)),
              mean_control_mu = mean(control_mu),
              mean_control_var = sd(control_mu) / sqrt(length(treat_mu))) %>%
    mutate(mean_treat_max = mean_treat_mu + mean_treat_var,
           mean_treat_min = mean_treat_mu - mean_treat_var,
           mean_control_max = mean_control_mu + mean_control_var,
           mean_control_min = mean_control_mu - mean_control_var)

df_proximal <- df_proximal %>%
    mutate(treat_mu = rowMeans(select(.,
                                      Control_1_S1_R1_001.sorted.bam,
                                      Control_2_S2_R1_001.sorted.bam,
                                      Control_3_S3_R1_001.sorted.bam)))
df_proximal <- df_proximal %>%
    mutate(control_mu = rowMeans(select(.,
                                        Taf_1_S4_R1_001.sorted.bam,
                                        Taf_2_S5_R1_001.sorted.bam,
                                        Taf_3_S6_R1_001.sorted.bam)))
df_proximal <- df_proximal %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                            Control_1_S1_R1_001.sorted.bam,
                                            Control_2_S2_R1_001.sorted.bam,
                                            Control_3_S3_R1_001.sorted.bam,
                                            Taf_1_S4_R1_001.sorted.bam,
                                            Taf_2_S5_R1_001.sorted.bam,
                                            Taf_3_S6_R1_001.sorted.bam))
final_proximal <- df_proximal %>% group_by(coord) %>%
    summarise(mean_treat_mu = mean(treat_mu),
              mean_treat_var = sd(treat_mu) / sqrt(length(treat_mu)),
              mean_control_mu = mean(control_mu),
              mean_control_var = sd(control_mu) / sqrt(length(treat_mu))) %>%
    mutate(mean_treat_max = mean_treat_mu + mean_treat_var,
           mean_treat_min = mean_treat_mu - mean_treat_var,
           mean_control_max = mean_control_mu + mean_control_var,
           mean_control_min = mean_control_mu - mean_control_var)


df_distal <- df_distal %>%
    mutate(treat_mu = rowMeans(select(.,
                                      Control_1_S1_R1_001.sorted.bam,
                                      Control_2_S2_R1_001.sorted.bam,
                                      Control_3_S3_R1_001.sorted.bam)))
df_distal <- df_distal %>%
    mutate(control_mu = rowMeans(select(.,
                                        Taf_1_S4_R1_001.sorted.bam,
                                        Taf_2_S5_R1_001.sorted.bam,
                                        Taf_3_S6_R1_001.sorted.bam)))
df_distal <- df_distal %>% subset(select = -c(geneid, Chr, Start, End, Length,
                                                    Control_1_S1_R1_001.sorted.bam,
                                                    Control_2_S2_R1_001.sorted.bam,
                                                    Control_3_S3_R1_001.sorted.bam,
                                                    Taf_1_S4_R1_001.sorted.bam,
                                                    Taf_2_S5_R1_001.sorted.bam,
                                                    Taf_3_S6_R1_001.sorted.bam))

final_distal <- df_distal %>% group_by(coord) %>%
    summarise(mean_treat_mu = mean(treat_mu),
              mean_treat_var = sd(treat_mu) / sqrt(length(treat_mu)),
              mean_control_mu = mean(control_mu),
              mean_control_var = sd(control_mu) / sqrt(length(treat_mu))) %>%
    mutate(mean_treat_max = mean_treat_mu + mean_treat_var,
           mean_treat_min = mean_treat_mu - mean_treat_var,
           mean_control_max = mean_control_mu + mean_control_var,
           mean_control_min = mean_control_mu - mean_control_var)

## Find the maximal coordinate for each
paused_treat_maxbin <- final_paused[which.max(final_paused$mean_treat_mu),]$coord
prox_treat_maxbin <- final_proximal[which.max(final_proximal$mean_treat_mu),]$coord
dist_treat_maxbin <- final_distal[which.max(final_distal$mean_treat_mu),]$coord

paused_control_maxbin <- final_paused[which.max(final_paused$mean_control_mu),]$coord
prox_control_maxbin <- final_proximal[which.max(final_proximal$mean_control_mu),]$coord
dist_control_maxbin <- final_distal[which.max(final_distal$mean_control_mu),]$coord

## Full Plot
library("ggthemes")
ggplot() + theme_tufte() +
    ## scale_color_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    ## scale_fill_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    ## Proximal Knockdown
    geom_line(data = final_proximal, aes(x = coord,
                                         y = mean_treat_mu,
                                      color = 'Proximal'), size = 0.25) +
    geom_ribbon(data = final_proximal, aes(x = coord,
                                           ymin = mean_treat_min,
                                        ymax = mean_treat_max,
                                        fill = 'Proximal'), size = 0.25, alpha = 0.2) +
    geom_vline(xintercept = prox_treat_maxbin, color = "#00BFC4", linetype = 'dashed') +
    ## Distal Knockdown
    geom_line(data = final_distal, aes(x = coord,
                                       y = mean_treat_mu,
                                       color = 'Distal'), size=0.25) +
    geom_ribbon(data = final_distal, aes(x = coord,
                                         ymin = mean_treat_min,
                                         ymax = mean_treat_max,
                                         fill = 'Distal'), size = 0.25, alpha = 0.2) +
    geom_vline(xintercept = dist_treat_maxbin, color = '#F8766D', linetype = 'dashed') +
    ## All Paused Genes
    geom_line(data = final_paused, aes(x = coord,
                                       y = mean_treat_mu,
                                       color = 'Paused'), size=0.25) +
    geom_ribbon(data = final_paused, aes(x = coord,
                                         ymin = mean_treat_min,
                                         ymax = mean_treat_max,
                                         fill = 'Paused'), size = 0.25, alpha = 0.2) +
    ## All Top Genes
    ## geom_line(data = final_top, aes(x = coord,
    ##                                 y = mean_treat_mu,
    ##                                 color = 'Top'), size=0.25) +
    ## geom_ribbon(data = final_top, aes(x = coord,
    ##                                   ymin = mean_treat_min,
    ##                                   ymax = mean_treat_max,
    ##                                   fill = 'Top'), size = 0.25, alpha = 0.2) +
    geom_hline(yintercept = 0, color = "black", size = 0.125) +
    labs(title = "Proximal vs Distal (Knockdown)", color = "Condition", fill = "Std. Dev. Mean") +
    xlab("Distance from TSS") + ylab("Normalized Read Depth") +
    ggsave("/home/zach/dowell_lab/pausing_meta_analysis/out/drosophila/metaplot/metaplot_knockdown.png", width = 10, height = 5)

ggplot() + theme_tufte() +
    ## scale_color_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    ## scale_fill_manual(values=c('Control'='#00BFC4', 'Knockdown'='#F8766D')) +
    ## Proximal Control
    geom_line(data = final_proximal, aes(x = coord,
                                         y = mean_control_mu,
                                      color = 'Proximal'), size = 0.25) +
    geom_ribbon(data = final_proximal, aes(x = coord,
                                           ymin = mean_control_min,
                                        ymax = mean_control_max,
                                        fill = 'Proximal'), size=0.25, alpha = 0.2) +
    geom_vline(xintercept = prox_control_maxbin, color = '#00BFC4', linetype = 'dashed') +
    ## Distal Control
    geom_line(data = final_distal, aes(x = coord,
                                       y = mean_control_mu,
                                       color = 'Distal'), size=0.25) +
    geom_ribbon(data = final_distal, aes(x = coord,
                                         ymin = mean_control_min,
                                            ymax = mean_control_max,
                                            fill = 'Distal'), size=0.25, alpha = 0.2) +
    ## All Paused Genes
    geom_line(data = final_paused, aes(x = coord,
                                       y = mean_control_mu,
                                       color = 'Paused'), size=0.25) +
    geom_ribbon(data = final_paused, aes(x = coord,
                                         ymin = mean_control_min,
                                         ymax = mean_control_max,
                                         fill = 'Paused'), size = 0.25, alpha = 0.2) +
    ## All Top Genes
    ## geom_line(data = final_top, aes(x = coord,
    ##                                 y = mean_control_mu,
    ##                                 color = 'Top'), size=0.25) +
    ## geom_ribbon(data = final_top, aes(x = coord,
    ##                                   ymin = mean_control_min,
    ##                                   ymax = mean_control_max,
    ##                                   fill = 'Top'), size = 0.25, alpha = 0.2) +
    geom_vline(xintercept = dist_control_maxbin, color = '#F8766D', linetype = 'dashed') +
    geom_hline(yintercept = 0, color = "black", size = 0.125) +
    labs(title = "Proximal vs Distal (Control)", color = "Condition", fill = "Std. Dev. Mean") +
    xlab("Distance from TSS") + ylab("Normalized Read Depth") +
    ggsave("/home/zach/dowell_lab/pausing_meta_analysis/out/drosophila/metaplot/metaplot_control.png", width = 10, height = 5)

## ggsave(outfile, width = 10, height = 5)

######################################################################
### metagene_graph_custom.r ends here
