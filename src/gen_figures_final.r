## This file contains code to generate preliminary figures for publication in our
## TAF1 PRO-Seq experiment, using the ggthemes and ggsci library to improve the
## appearance of those figures.  Additionally, we filter for a set of well
## expressed genes as defined by Ebmeier 2017
suppressMessages(library("tidyverse"))
library(ggthemes)
library(ggsci)

## Well expressed genes
wellexpr <- read_delim("expressed_genes_filtertable.txt", col_names = c("name", "tx_name"), 
    delim = "\t")

## Generate figures...
c_1_PROD <- read_delim("C413_1_S3_R1_001.trim.bedGraph_pause_ratios_PROD.data", col_names = c("tx_name", 
    "strand", "c_1_PROD", "coverage_c1_PROD"), delim = " ") %>% subset(select = -strand)
c_2_PROD <- read_delim("C413_2_S4_R1_001.trim.bedGraph_pause_ratios_PROD.data", col_names = c("tx_name", 
    "strand", "c_2_PROD", "coverage_c2_PROD"), delim = " ") %>% subset(select = -strand)
p_1_PROD <- read_delim("PO_1_S1_R1_001.trim.bedGraph_pause_ratios_PROD.data", col_names = c("tx_name", 
    "strand", "p_1_PROD", "coverage_p1_PROD"), delim = " ") %>% subset(select = -strand)
p_2_PROD <- read_delim("PO_2_S2_R1_001.trim.bedGraph_pause_ratios_PROD.data", col_names = c("tx_name", 
    "strand", "p_2_PROD", "coverage_p2_PROD"), delim = " ") %>% subset(select = -strand)

## Calculate statistics for each gene in our list
ddt_PROD <- list(c_1_PROD, c_2_PROD, p_1_PROD, p_2_PROD) %>% reduce(left_join, by = "tx_name")
ddt_PROD$c_pause_mean <- rowMeans(ddt_PROD[c("c_1_PROD", "c_2_PROD")], na.rm = TRUE)
ddt_PROD$p_pause_mean <- rowMeans(ddt_PROD[c("p_1_PROD", "p_2_PROD")], na.rm = TRUE)
ddt_PROD$c_coverage_mean <- rowMeans(ddt_PROD[c("coverage_c1_PROD", "coverage_c2_PROD")], 
    na.rm = FALSE)
ddt_PROD$p_coverage_mean <- rowMeans(ddt_PROD[c("coverage_p1_PROD", "coverage_p2_PROD")], 
    na.rm = FALSE)

## ddt_PROD <- subset(ddt_PROD, !duplicated(tx_name))

## Calculate the log2 difference between the means of the two samples
ddt_PROD <- ddt_PROD %>% mutate(pause_diff = log2(c_pause_mean) - log2(p_pause_mean))
ddt_PROD$coverage_mean <- rowMeans(ddt_PROD[c("c_coverage_mean", "p_coverage_mean")], 
    na.rm = TRUE)

## Filter out well expressed genes according to our filter table.
ddt_PROD_wellexpressed <- subset(ddt_PROD, tx_name %in% wellexpr$tx_name)

## scatter plot
ggplot(data = ddt_PROD, mapping = aes(x = p_pause_mean, y = c_pause_mean), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + 
    scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/prepublish_scatter.png"), plot = last_plot(), 
    device = "png")

## ecdf plot
ggplot(data = ddt_PROD) + stat_ecdf(aes(p_pause_mean, color = "Control")) + stat_ecdf(aes(c_pause_mean, 
    color = "Treatment")) + theme_tufte() + scale_color_npg() + theme(legend.title = element_blank()) + 
    labs(x = "Pausing Index", y = "Cumulative Distribution", title = "Cumulative Distribution") + 
    scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/prepublish_ecdf.png"), plot = last_plot(), 
    device = "png")

fddt_PROD <- filter(ddt_PROD, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_PROD, mapping = aes(x = coverage_mean, y = pause_diff), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + theme_tufte() + 
    labs(x = "Normalized Coverage Level", y = "Change in Pause Index", title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + 
    scale_x_log10() + ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/prepublish_MA.png"), plot = last_plot(), 
    device = "png")

############### WELL EXPRESSED ################################

## PROD_wellexpressed scatter
ggplot(data = ddt_PROD_wellexpressed, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/prepublish_wellexpressed_scatter.png"), 
    plot = last_plot(), device = "png")

## PROD_wellexpressed ecdf
ggplot(data = ddt_PROD_wellexpressed) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/prepublish_wellexpressed_ecdf.png"), 
    plot = last_plot(), device = "png")

fddt_PROD_wellexpressed <- filter(ddt_PROD_wellexpressed, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_PROD_wellexpressed, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/prepublish_wellexpressed_MA.png"), 
    plot = last_plot(), device = "png")
