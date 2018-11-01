suppressMessages(library("tidyverse"))
library(ggthemes)
library(ggsci)

## Well expressed genes
wellexpr <- read_delim("expressed_genes_filtertable.txt", col_names=c("name", 'tx_name'), delim="\t")
inr <- read_delim("hg38_matched_genes_BBCABW.data", col_names=c("name", 'sequence'), delim="\t")
mte <- read_delim("hg38_matched_genes_CGANC....CGG.data", col_names=c("name", 'sequence'), delim="\t")
dpe <- read_delim("hg38_matched_genes_RGWYVT.data", col_names=c("name", 'sequence'), delim="\t")
tata <- read_delim("hg38_matched_genes_WWWW.data", col_names=c("name", 'sequence'), delim="\t")

## Generate figures...
c_1_301 <- read_delim("C413_1_S3_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                     col_names=c('tx_name', 'strand', 'c_1_301', 'coverage_c1_301'), delim=" ") %>%
    subset(select=-strand)
c_2_301 <- read_delim("C413_2_S4_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'c_2_301', 'coverage_c2_301'), delim=" ") %>%
    subset(select=-strand)
p_1_301 <- read_delim("PO_1_S1_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'p_1_301', 'coverage_p1_301'), delim=" ") %>%
    subset(select=-strand)
p_2_301 <- read_delim("PO_2_S2_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'p_2_301', 'coverage_p2_301'), delim=" ") %>%
    subset(select=-strand)

ddt_301 <- list(c_1_301, c_2_301, p_1_301, p_2_301) %>% reduce(left_join, by="tx_name")
ddt_301$c_pause_mean <- rowMeans(ddt_301[c('c_1_301', 'c_2_301')], na.rm=TRUE)
ddt_301$p_pause_mean <- rowMeans(ddt_301[c('p_1_301', 'p_2_301')], na.rm=TRUE)
ddt_301$c_coverage_mean <- rowMeans(ddt_301[c('coverage_c1_301', 'coverage_c2_301')], na.rm=FALSE)
ddt_301$p_coverage_mean <- rowMeans(ddt_301[c('coverage_p1_301', 'coverage_p2_301')], na.rm=FALSE)

ddt_301 <- subset(ddt_301, !duplicated(tx_name))
ddt_301 <- ddt_301 %>% mutate(pause_diff = log2(c_pause_mean) - log2(p_pause_mean))
ddt_301$coverage_mean <- rowMeans(ddt_301[c('c_coverage_mean', 'p_coverage_mean')], na.rm=TRUE)

ddt_301_wellexpressed <- subset(ddt_301, tx_name %in% wellexpr$tx_name)
ddt_301_wellexpressed_inr <- subset(ddt_301_wellexpressed, tx_name %in% inr$name)
ddt_301_wellexpressed_mte <- subset(ddt_301_wellexpressed, tx_name %in% mte$name)
ddt_301_wellexpressed_dpe <- subset(ddt_301_wellexpressed, tx_name %in% dpe$name)
ddt_301_wellexpressed_tata <- subset(ddt_301_wellexpressed, tx_name %in% tata$name)

ddt_301_inr <- subset(ddt_301, tx_name %in% inr$name)
ddt_301_mte <- subset(ddt_301, tx_name %in% mte$name)
ddt_301_dpe <- subset(ddt_301, tx_name %in% dpe$name)
ddt_301_tata <- subset(ddt_301, tx_name %in% tata$name)

ddt_301_inr_and_mte <- subset(ddt_301, tx_name %in% inr$name & tx_name %in% mte$name)
ddt_301_inr_and_dpe <- subset(ddt_301, tx_name %in% inr$name & tx_name %in% dpe$name)
ddt_301_inr_and_tata <- subset(ddt_301, tx_name %in% inr$name & tx_name %in% tata$name)
ddt_301_tata_and_mte <- subset(ddt_301, tx_name %in% tata$name & tx_name %in% mte$name)
ddt_301_tata_and_dpe <- subset(ddt_301, tx_name %in% tata$name & tx_name %in% dpe$name)

ddt_301_wellexpressed_inr_and_mte <- subset(ddt_301_wellexpressed, tx_name %in% inr$name & tx_name %in% mte$name)
ddt_301_wellexpressed_inr_and_dpe <- subset(ddt_301_wellexpressed, tx_name %in% inr$name & tx_name %in% dpe$name)
ddt_301_wellexpressed_inr_and_tata <- subset(ddt_301_wellexpressed, tx_name %in% inr$name & tx_name %in% tata$name)
ddt_301_wellexpressed_tata_and_mte <- subset(ddt_301_wellexpressed, tx_name %in% tata$name & tx_name %in% mte$name)
ddt_301_wellexpressed_tata_and_dpe <- subset(ddt_301_wellexpressed, tx_name %in% tata$name & tx_name %in% dpe$name)

ddt_301_all <- subset(ddt_301, tx_name %in% inr$name & tx_name %in% tata$name &
                              tx_name %in% dpe$name & tx_name %in% mte$name)

## 301 scatter
ggplot(data = ddt_301, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
theme_tufte() +
scale_fill_gsea() +
labs(x = "Control",
     y= "Treatment",
     title= "Change in Pausing Index") +
scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_scatter.png"),
    plot = last_plot(), device = "png")

## 301 ecdf
ggplot(data = ddt_301) +
stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
theme_tufte() +
scale_color_npg() +
theme(legend.title=element_blank()) +
labs(x = "Pausing Index",
     y= "Cumulative Distribution",
     title= "Cumulative Distribution") +
scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301 <- filter(ddt_301, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
geom_bin2d(bins=75) + geom_hline(yintercept=0) +
scale_fill_gsea() +
theme_tufte() +
labs(x = "Normalized Coverage Level",
     y= "Change in Pause Index",
     title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
scale_x_log10() +
## scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_MA.png"),
    plot = last_plot(), device = "png")

################################################################################
################################################################################
################################################################################

## 301_wellexpressed scatter
ggplot(data = ddt_301_wellexpressed, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
theme_tufte() +
scale_fill_gsea() +
labs(x = "Control",
     y= "Treatment",
     title= "Change in Pausing Index") +
scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed ecdf
ggplot(data = ddt_301_wellexpressed) +
stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
theme_tufte() +
scale_color_npg() +
theme(legend.title=element_blank()) +
labs(x = "Pausing Index",
     y= "Cumulative Distribution",
     title= "Cumulative Distribution") +
scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_wellexpressed <- filter(ddt_301_wellexpressed, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_wellexpressed, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
geom_bin2d(bins=75) + geom_hline(yintercept=0) +
scale_fill_gsea() +
theme_tufte() +
labs(x = "Normalized Coverage Level",
     y= "Change in Pause Index",
     title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
scale_x_log10() +
## scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_MA.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_tata scatter
ggplot(data = ddt_301_wellexpressed_tata, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
theme_tufte() +
scale_fill_gsea() +
labs(x = "Control",
     y= "Treatment",
     title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_tata ecdf
ggplot(data = ddt_301_wellexpressed_tata) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_wellexpressed_tata <- filter(ddt_301_wellexpressed_tata, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_tata, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_MA.png"),
        plot = last_plot(), device = "png")

## 301_wellexpressed_dpe scatter
ggplot(data = ddt_301_wellexpressed_dpe, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_dpe_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_dpe ecdf
ggplot(data = ddt_301_wellexpressed_dpe) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_dpe_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_wellexpressed_dpe <- filter(ddt_301_wellexpressed_dpe, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_dpe, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_dpe_MA.png"),
        plot = last_plot(), device = "png")

## 301_wellexpressed_mte scatter
ggplot(data = ddt_301_wellexpressed_mte, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_mte_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_mte ecdf
ggplot(data = ddt_301_wellexpressed_mte) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_mte_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_wellexpressed_mte <- filter(ddt_301_wellexpressed_mte, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_mte, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_mte_MA.png"),
        plot = last_plot(), device = "png")

## 301_wellexpressed_inr scatter
ggplot(data = ddt_301_wellexpressed_inr, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_inr ecdf
ggplot(data = ddt_301_wellexpressed_inr) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_wellexpressed_inr <- filter(ddt_301_wellexpressed_inr, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_inr, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_MA.png"),
        plot = last_plot(), device = "png")

################################################################################
################################################################################
################################################################################

## 301_tata scatter
ggplot(data = ddt_301_tata, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_tata_scatter.png"),
    plot = last_plot(), device = "png")

## 301_tata ecdf
ggplot(data = ddt_301_tata) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_tata_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_tata <- filter(ddt_301_tata, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_tata, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_tata_MA.png"),
        plot = last_plot(), device = "png")

## 301_dpe scatter
ggplot(data = ddt_301_dpe, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_dpe_scatter.png"),
    plot = last_plot(), device = "png")

## 301_dpe ecdf
ggplot(data = ddt_301_dpe) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_dpe_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_dpe <- filter(ddt_301_dpe, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_dpe, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_dpe_MA.png"),
        plot = last_plot(), device = "png")

## 301_mte scatter
ggplot(data = ddt_301_mte, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_mte_scatter.png"),
    plot = last_plot(), device = "png")

## 301_mte ecdf
ggplot(data = ddt_301_mte) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_mte_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_mte <- filter(ddt_301_mte, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_mte, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_mte_MA.png"),
        plot = last_plot(), device = "png")

## 301_inr scatter
ggplot(data = ddt_301_inr, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_scatter.png"),
    plot = last_plot(), device = "png")

## 301_inr ecdf
ggplot(data = ddt_301_inr) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301_inr <- filter(ddt_301_inr, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301_inr, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_inr_MA.png"),
        plot = last_plot(), device = "png")

################################################################################
################################################################################
################################################################################

## 301_inr_and_mte scatter
ggplot(data = ddt_301_inr_and_mte, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_and_mte_scatter.png"),
    plot = last_plot(), device = "png")

## 301_inr_and_mte ecdf
ggplot(data = ddt_301_inr_and_mte) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_and_mte_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_inr_and_mte, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_inr_and_mte_MA.png"),
        plot = last_plot(), device = "png")

## 301_inr_and_dpe scatter
ggplot(data = ddt_301_inr_and_dpe, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_and_dpe_scatter.png"),
    plot = last_plot(), device = "png")

## 301_inr_and_dpe ecdf
ggplot(data = ddt_301_inr_and_dpe) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_and_dpe_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_inr_and_dpe, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_inr_and_dpe_MA.png"),
        plot = last_plot(), device = "png")

## 301_inr_and_tata scatter
ggplot(data = ddt_301_inr_and_tata, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_and_tata_scatter.png"),
    plot = last_plot(), device = "png")

## 301_inr_and_tata ecdf
ggplot(data = ddt_301_inr_and_tata) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_inr_and_tata_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_inr_and_tata, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_inr_and_tata_MA.png"),
        plot = last_plot(), device = "png")

## 301_tata_and_mte scatter
ggplot(data = ddt_301_tata_and_mte, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_tata_and_mte_scatter.png"),
    plot = last_plot(), device = "png")

## 301_tata_and_mte ecdf
ggplot(data = ddt_301_tata_and_mte) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_tata_and_mte_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_tata_and_mte, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_tata_and_mte_MA.png"),
        plot = last_plot(), device = "png")

## 301_tata_and_dpe scatter
ggplot(data = ddt_301_tata_and_dpe, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_tata_and_dpe_scatter.png"),
    plot = last_plot(), device = "png")

## 301_tata_and_dpe ecdf
ggplot(data = ddt_301_tata_and_dpe) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_tata_and_dpe_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_tata_and_dpe, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_tata_and_dpe_MA.png"),
        plot = last_plot(), device = "png")

################################################################################
################################################################################
################################################################################

## 301_wellexpressed_inr_and_mte scatter
ggplot(data = ddt_301_wellexpressed_inr_and_mte, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_mte_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_inr_and_mte ecdf
ggplot(data = ddt_301_wellexpressed_inr_and_mte) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_mte_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_inr_and_mte, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_mte_MA.png"),
        plot = last_plot(), device = "png")

## 301_wellexpressed_inr_and_dpe scatter
ggplot(data = ddt_301_wellexpressed_inr_and_dpe, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_dpe_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_inr_and_dpe ecdf
ggplot(data = ddt_301_wellexpressed_inr_and_dpe) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_dpe_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_inr_and_dpe, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_dpe_MA.png"),
        plot = last_plot(), device = "png")

## 301_wellexpressed_inr_and_tata scatter
ggplot(data = ddt_301_wellexpressed_inr_and_tata, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_tata_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_inr_and_tata ecdf
ggplot(data = ddt_301_wellexpressed_inr_and_tata) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_tata_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_inr_and_tata, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_inr_and_tata_MA.png"),
        plot = last_plot(), device = "png")

## 301_wellexpressed_tata_and_mte scatter
ggplot(data = ddt_301_wellexpressed_tata_and_mte, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_and_mte_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_tata_and_mte ecdf
ggplot(data = ddt_301_wellexpressed_tata_and_mte) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_and_mte_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_tata_and_mte, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_and_mte_MA.png"),
        plot = last_plot(), device = "png")

## 301_wellexpressed_tata_and_dpe scatter
ggplot(data = ddt_301_wellexpressed_tata_and_dpe, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
    geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    scale_fill_gsea() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_and_dpe_scatter.png"),
    plot = last_plot(), device = "png")

## 301_wellexpressed_tata_and_dpe ecdf
ggplot(data = ddt_301_wellexpressed_tata_and_dpe) +
    stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
    theme_tufte() +
    scale_color_npg() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_and_dpe_ecdf.png"),
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_301_wellexpressed_tata_and_dpe, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
    geom_bin2d(bins=75) + geom_hline(yintercept=0) +
    scale_fill_gsea() +
    theme_tufte() +
    labs(x = "Normalized Coverage Level",
         y= "Change in Pause Index",
         title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
    scale_x_log10() +
    ## scale_y_log10()
    ggsave(
        str_c("/scratch/Users/zama8258/pause_output/301_wellexpressed_tata_and_dpe_MA.png"),
        plot = last_plot(), device = "png")

################################################################################
################################################################################
################################################################################
