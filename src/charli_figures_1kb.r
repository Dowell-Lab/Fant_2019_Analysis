suppressMessages(library("tidyverse"))
library(ggthemes)
library(ggsci)

## Well expressed genes
wellexpr <- read_delim("expressed_genes_filtertable.txt", col_names = c("name", "tx_name"), 
    delim = "\t")
inr <- read_delim("hg38_matched_genes_BBCABW.data", col_names = c("name", "sequence"), 
    delim = "\t")
mte <- read_delim("hg38_matched_genes_CGANC....CGG.data", col_names = c("name", "sequence"), 
    delim = "\t")
dpe <- read_delim("hg38_matched_genes_RGWYVT.data", col_names = c("name", "sequence"), 
    delim = "\t")
tata <- read_delim("hg38_matched_genes_WWWW.data", col_names = c("name", "sequence"), 
    delim = "\t")

## Generate figures...
c_1_1000 <- read_delim("C413_1_S3_R1_001.trim.bedGraph_pause_ratios_1000.data", col_names = c("tx_name", 
    "strand", "c_1_1000", "coverage_c1_1000"), delim = " ") %>% subset(select = -strand)
c_2_1000 <- read_delim("C413_2_S4_R1_001.trim.bedGraph_pause_ratios_1000.data", col_names = c("tx_name", 
    "strand", "c_2_1000", "coverage_c2_1000"), delim = " ") %>% subset(select = -strand)
p_1_1000 <- read_delim("PO_1_S1_R1_001.trim.bedGraph_pause_ratios_1000.data", col_names = c("tx_name", 
    "strand", "p_1_1000", "coverage_p1_1000"), delim = " ") %>% subset(select = -strand)
p_2_1000 <- read_delim("PO_2_S2_R1_001.trim.bedGraph_pause_ratios_1000.data", col_names = c("tx_name", 
    "strand", "p_2_1000", "coverage_p2_1000"), delim = " ") %>% subset(select = -strand)

ddt_1000 <- list(c_1_1000, c_2_1000, p_1_1000, p_2_1000) %>% reduce(left_join, by = "tx_name")
ddt_1000$c_pause_mean <- rowMeans(ddt_1000[c("c_1_1000", "c_2_1000")], na.rm = TRUE)
ddt_1000$p_pause_mean <- rowMeans(ddt_1000[c("p_1_1000", "p_2_1000")], na.rm = TRUE)
ddt_1000$c_coverage_mean <- rowMeans(ddt_1000[c("coverage_c1_1000", "coverage_c2_1000")], 
    na.rm = FALSE)
ddt_1000$p_coverage_mean <- rowMeans(ddt_1000[c("coverage_p1_1000", "coverage_p2_1000")], 
    na.rm = FALSE)

ddt_1000 <- subset(ddt_1000, !duplicated(tx_name))
ddt_1000 <- ddt_1000 %>% mutate(pause_diff = log2(c_pause_mean) - log2(p_pause_mean))
ddt_1000$coverage_mean <- rowMeans(ddt_1000[c("c_coverage_mean", "p_coverage_mean")], 
    na.rm = TRUE)

ddt_1000_wellexpressed <- subset(ddt_1000, tx_name %in% wellexpr$tx_name)
ddt_1000_wellexpressed_inr <- subset(ddt_1000_wellexpressed, tx_name %in% inr$name)
ddt_1000_wellexpressed_mte <- subset(ddt_1000_wellexpressed, tx_name %in% mte$name)
ddt_1000_wellexpressed_dpe <- subset(ddt_1000_wellexpressed, tx_name %in% dpe$name)
ddt_1000_wellexpressed_tata <- subset(ddt_1000_wellexpressed, tx_name %in% tata$name)

ddt_1000_inr <- subset(ddt_1000, tx_name %in% inr$name)
ddt_1000_mte <- subset(ddt_1000, tx_name %in% mte$name)
ddt_1000_dpe <- subset(ddt_1000, tx_name %in% dpe$name)
ddt_1000_tata <- subset(ddt_1000, tx_name %in% tata$name)

ddt_1000_inr_and_mte <- subset(ddt_1000, tx_name %in% inr$name & tx_name %in% mte$name)
ddt_1000_inr_and_dpe <- subset(ddt_1000, tx_name %in% inr$name & tx_name %in% dpe$name)
ddt_1000_inr_and_tata <- subset(ddt_1000, tx_name %in% inr$name & tx_name %in% tata$name)
ddt_1000_tata_and_mte <- subset(ddt_1000, tx_name %in% tata$name & tx_name %in% mte$name)
ddt_1000_tata_and_dpe <- subset(ddt_1000, tx_name %in% tata$name & tx_name %in% dpe$name)

ddt_1000_wellexpressed_inr_and_mte <- subset(ddt_1000_wellexpressed, tx_name %in% 
    inr$name & tx_name %in% mte$name)
ddt_1000_wellexpressed_inr_and_dpe <- subset(ddt_1000_wellexpressed, tx_name %in% 
    inr$name & tx_name %in% dpe$name)
ddt_1000_wellexpressed_inr_and_tata <- subset(ddt_1000_wellexpressed, tx_name %in% 
    inr$name & tx_name %in% tata$name)
ddt_1000_wellexpressed_tata_and_mte <- subset(ddt_1000_wellexpressed, tx_name %in% 
    tata$name & tx_name %in% mte$name)
ddt_1000_wellexpressed_tata_and_dpe <- subset(ddt_1000_wellexpressed, tx_name %in% 
    tata$name & tx_name %in% dpe$name)

ddt_1000_all <- subset(ddt_1000, tx_name %in% inr$name & tx_name %in% tata$name & 
    tx_name %in% dpe$name & tx_name %in% mte$name)

## 1000 scatter
ggplot(data = ddt_1000, mapping = aes(x = p_pause_mean, y = c_pause_mean), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + 
    scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_scatter.png"), plot = last_plot(), 
    device = "png")

## 1000 ecdf
ggplot(data = ddt_1000) + stat_ecdf(aes(p_pause_mean, color = "Control")) + stat_ecdf(aes(c_pause_mean, 
    color = "Treatment")) + theme_tufte() + scale_color_npg() + theme(legend.title = element_blank()) + 
    labs(x = "Pausing Index", y = "Cumulative Distribution", title = "Cumulative Distribution") + 
    scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_ecdf.png"), plot = last_plot(), 
    device = "png")

fddt_1000 <- filter(ddt_1000, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000, mapping = aes(x = coverage_mean, y = pause_diff), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + theme_tufte() + 
    labs(x = "Normalized Coverage Level", y = "Change in Pause Index", title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + 
    scale_x_log10() + ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_MA.png"), plot = last_plot(), 
    device = "png")

################################################################################ 

## 1000_wellexpressed scatter
ggplot(data = ddt_1000_wellexpressed, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed ecdf
ggplot(data = ddt_1000_wellexpressed) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_ecdf.png"), 
    plot = last_plot(), device = "png")

fddt_1000_wellexpressed <- filter(ddt_1000_wellexpressed, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_wellexpressed_tata scatter
ggplot(data = ddt_1000_wellexpressed_tata, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_tata ecdf
ggplot(data = ddt_1000_wellexpressed_tata) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_ecdf.png"), 
    plot = last_plot(), device = "png")

fddt_1000_wellexpressed_tata <- filter(ddt_1000_wellexpressed_tata, pause_diff > 
    -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_tata, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_MA.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_dpe scatter
ggplot(data = ddt_1000_wellexpressed_dpe, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_dpe_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_dpe ecdf
ggplot(data = ddt_1000_wellexpressed_dpe) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_dpe_ecdf.png"), 
    plot = last_plot(), device = "png")

fddt_1000_wellexpressed_dpe <- filter(ddt_1000_wellexpressed_dpe, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_dpe, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_dpe_MA.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_mte scatter
ggplot(data = ddt_1000_wellexpressed_mte, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_mte_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_mte ecdf
ggplot(data = ddt_1000_wellexpressed_mte) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_mte_ecdf.png"), 
    plot = last_plot(), device = "png")

fddt_1000_wellexpressed_mte <- filter(ddt_1000_wellexpressed_mte, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_mte, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_mte_MA.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_inr scatter
ggplot(data = ddt_1000_wellexpressed_inr, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_inr ecdf
ggplot(data = ddt_1000_wellexpressed_inr) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_ecdf.png"), 
    plot = last_plot(), device = "png")

fddt_1000_wellexpressed_inr <- filter(ddt_1000_wellexpressed_inr, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_inr, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_MA.png"), 
    plot = last_plot(), device = "png")

################################################################################ 

## 1000_tata scatter
ggplot(data = ddt_1000_tata, mapping = aes(x = p_pause_mean, y = c_pause_mean), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + 
    scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_scatter.png"), plot = last_plot(), 
    device = "png")

## 1000_tata ecdf
ggplot(data = ddt_1000_tata) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_ecdf.png"), plot = last_plot(), 
    device = "png")

fddt_1000_tata <- filter(ddt_1000_tata, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_tata, mapping = aes(x = coverage_mean, y = pause_diff), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + theme_tufte() + 
    labs(x = "Normalized Coverage Level", y = "Change in Pause Index", title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + 
    scale_x_log10() + ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_dpe scatter
ggplot(data = ddt_1000_dpe, mapping = aes(x = p_pause_mean, y = c_pause_mean), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + 
    scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_dpe_scatter.png"), plot = last_plot(), 
    device = "png")

## 1000_dpe ecdf
ggplot(data = ddt_1000_dpe) + stat_ecdf(aes(p_pause_mean, color = "Control")) + stat_ecdf(aes(c_pause_mean, 
    color = "Treatment")) + theme_tufte() + scale_color_npg() + theme(legend.title = element_blank()) + 
    labs(x = "Pausing Index", y = "Cumulative Distribution", title = "Cumulative Distribution") + 
    scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_dpe_ecdf.png"), plot = last_plot(), 
    device = "png")

fddt_1000_dpe <- filter(ddt_1000_dpe, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_dpe, mapping = aes(x = coverage_mean, y = pause_diff), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + theme_tufte() + 
    labs(x = "Normalized Coverage Level", y = "Change in Pause Index", title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + 
    scale_x_log10() + ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_dpe_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_mte scatter
ggplot(data = ddt_1000_mte, mapping = aes(x = p_pause_mean, y = c_pause_mean), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + 
    scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_mte_scatter.png"), plot = last_plot(), 
    device = "png")

## 1000_mte ecdf
ggplot(data = ddt_1000_mte) + stat_ecdf(aes(p_pause_mean, color = "Control")) + stat_ecdf(aes(c_pause_mean, 
    color = "Treatment")) + theme_tufte() + scale_color_npg() + theme(legend.title = element_blank()) + 
    labs(x = "Pausing Index", y = "Cumulative Distribution", title = "Cumulative Distribution") + 
    scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_mte_ecdf.png"), plot = last_plot(), 
    device = "png")

fddt_1000_mte <- filter(ddt_1000_mte, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_mte, mapping = aes(x = coverage_mean, y = pause_diff), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + theme_tufte() + 
    labs(x = "Normalized Coverage Level", y = "Change in Pause Index", title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + 
    scale_x_log10() + ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_mte_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_inr scatter
ggplot(data = ddt_1000_inr, mapping = aes(x = p_pause_mean, y = c_pause_mean), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + theme_tufte() + 
    scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_scatter.png"), plot = last_plot(), 
    device = "png")

## 1000_inr ecdf
ggplot(data = ddt_1000_inr) + stat_ecdf(aes(p_pause_mean, color = "Control")) + stat_ecdf(aes(c_pause_mean, 
    color = "Treatment")) + theme_tufte() + scale_color_npg() + theme(legend.title = element_blank()) + 
    labs(x = "Pausing Index", y = "Cumulative Distribution", title = "Cumulative Distribution") + 
    scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_ecdf.png"), plot = last_plot(), 
    device = "png")

fddt_1000_inr <- filter(ddt_1000_inr, pause_diff > -6e-05)

## MA-ish plot
ggplot(data = ddt_1000_inr, mapping = aes(x = coverage_mean, y = pause_diff), alpha = 1/10) + 
    geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + theme_tufte() + 
    labs(x = "Normalized Coverage Level", y = "Change in Pause Index", title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + 
    scale_x_log10() + ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_MA.png"), plot = last_plot(), 
    device = "png")

################################################################################ 

## 1000_inr_and_mte scatter
ggplot(data = ddt_1000_inr_and_mte, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_mte_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_inr_and_mte ecdf
ggplot(data = ddt_1000_inr_and_mte) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_mte_ecdf.png"), plot = last_plot(), 
    device = "png")

## MA-ish plot
ggplot(data = ddt_1000_inr_and_mte, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_mte_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_inr_and_dpe scatter
ggplot(data = ddt_1000_inr_and_dpe, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_dpe_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_inr_and_dpe ecdf
ggplot(data = ddt_1000_inr_and_dpe) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_dpe_ecdf.png"), plot = last_plot(), 
    device = "png")

## MA-ish plot
ggplot(data = ddt_1000_inr_and_dpe, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_dpe_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_inr_and_tata scatter
ggplot(data = ddt_1000_inr_and_tata, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_tata_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_inr_and_tata ecdf
ggplot(data = ddt_1000_inr_and_tata) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_tata_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_inr_and_tata, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_inr_and_tata_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_tata_and_mte scatter
ggplot(data = ddt_1000_tata_and_mte, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_and_mte_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_tata_and_mte ecdf
ggplot(data = ddt_1000_tata_and_mte) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_and_mte_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_tata_and_mte, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_and_mte_MA.png"), plot = last_plot(), 
    device = "png")

## 1000_tata_and_dpe scatter
ggplot(data = ddt_1000_tata_and_dpe, mapping = aes(x = p_pause_mean, y = c_pause_mean), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, slope = 1)) + 
    theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", title = "Change in Pausing Index") + 
    scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_and_dpe_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_tata_and_dpe ecdf
ggplot(data = ddt_1000_tata_and_dpe) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_and_dpe_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_tata_and_dpe, mapping = aes(x = coverage_mean, y = pause_diff), 
    alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + scale_fill_gsea() + 
    theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_tata_and_dpe_MA.png"), plot = last_plot(), 
    device = "png")

################################################################################ 

## 1000_wellexpressed_inr_and_mte scatter
ggplot(data = ddt_1000_wellexpressed_inr_and_mte, mapping = aes(x = p_pause_mean, 
    y = c_pause_mean), alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, 
    slope = 1)) + theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", 
    title = "Change in Pausing Index") + scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_mte_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_inr_and_mte ecdf
ggplot(data = ddt_1000_wellexpressed_inr_and_mte) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_mte_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_inr_and_mte, mapping = aes(x = coverage_mean, 
    y = pause_diff), alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + 
    scale_fill_gsea() + theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_mte_MA.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_inr_and_dpe scatter
ggplot(data = ddt_1000_wellexpressed_inr_and_dpe, mapping = aes(x = p_pause_mean, 
    y = c_pause_mean), alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, 
    slope = 1)) + theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", 
    title = "Change in Pausing Index") + scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_dpe_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_inr_and_dpe ecdf
ggplot(data = ddt_1000_wellexpressed_inr_and_dpe) + stat_ecdf(aes(p_pause_mean, color = "Control")) + 
    stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + scale_color_npg() + 
    theme(legend.title = element_blank()) + labs(x = "Pausing Index", y = "Cumulative Distribution", 
    title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_dpe_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_inr_and_dpe, mapping = aes(x = coverage_mean, 
    y = pause_diff), alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + 
    scale_fill_gsea() + theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_dpe_MA.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_inr_and_tata scatter
ggplot(data = ddt_1000_wellexpressed_inr_and_tata, mapping = aes(x = p_pause_mean, 
    y = c_pause_mean), alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, 
    slope = 1)) + theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", 
    title = "Change in Pausing Index") + scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_tata_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_inr_and_tata ecdf
ggplot(data = ddt_1000_wellexpressed_inr_and_tata) + stat_ecdf(aes(p_pause_mean, 
    color = "Control")) + stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + 
    scale_color_npg() + theme(legend.title = element_blank()) + labs(x = "Pausing Index", 
    y = "Cumulative Distribution", title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_tata_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_inr_and_tata, mapping = aes(x = coverage_mean, 
    y = pause_diff), alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + 
    scale_fill_gsea() + theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_inr_and_tata_MA.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_tata_and_mte scatter
ggplot(data = ddt_1000_wellexpressed_tata_and_mte, mapping = aes(x = p_pause_mean, 
    y = c_pause_mean), alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, 
    slope = 1)) + theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", 
    title = "Change in Pausing Index") + scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_and_mte_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_tata_and_mte ecdf
ggplot(data = ddt_1000_wellexpressed_tata_and_mte) + stat_ecdf(aes(p_pause_mean, 
    color = "Control")) + stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + 
    scale_color_npg() + theme(legend.title = element_blank()) + labs(x = "Pausing Index", 
    y = "Cumulative Distribution", title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_and_mte_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_tata_and_mte, mapping = aes(x = coverage_mean, 
    y = pause_diff), alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + 
    scale_fill_gsea() + theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_and_mte_MA.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_tata_and_dpe scatter
ggplot(data = ddt_1000_wellexpressed_tata_and_dpe, mapping = aes(x = p_pause_mean, 
    y = c_pause_mean), alpha = 1/10) + geom_bin2d(bins = 75) + geom_abline(aes(intercept = 0, 
    slope = 1)) + theme_tufte() + scale_fill_gsea() + labs(x = "Control", y = "Treatment", 
    title = "Change in Pausing Index") + scale_x_log10() + scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_and_dpe_scatter.png"), 
    plot = last_plot(), device = "png")

## 1000_wellexpressed_tata_and_dpe ecdf
ggplot(data = ddt_1000_wellexpressed_tata_and_dpe) + stat_ecdf(aes(p_pause_mean, 
    color = "Control")) + stat_ecdf(aes(c_pause_mean, color = "Treatment")) + theme_tufte() + 
    scale_color_npg() + theme(legend.title = element_blank()) + labs(x = "Pausing Index", 
    y = "Cumulative Distribution", title = "Cumulative Distribution") + scale_x_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_and_dpe_ecdf.png"), 
    plot = last_plot(), device = "png")

## MA-ish plot
ggplot(data = ddt_1000_wellexpressed_tata_and_dpe, mapping = aes(x = coverage_mean, 
    y = pause_diff), alpha = 1/10) + geom_bin2d(bins = 75) + geom_hline(yintercept = 0) + 
    scale_fill_gsea() + theme_tufte() + labs(x = "Normalized Coverage Level", y = "Change in Pause Index", 
    title = "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") + scale_x_log10() + 
    ## scale_y_log10()
ggsave(str_c("/scratch/Users/zama8258/pause_output/1000_wellexpressed_tata_and_dpe_MA.png"), 
    plot = last_plot(), device = "png")

################################################################################ 
