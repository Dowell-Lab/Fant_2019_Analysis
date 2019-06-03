suppressMessages(library("tidyverse"))
library(ggthemes)
library(ggsci)

# TODO proximal and distal gense from LIS

setwd("/home/zach/dowell_lab/pausing_meta_analysis/out/drosophila/pausing")
## Generate figures...
## You only need to adjust the filenames here...
c_1_301 <- read_delim("Control_1_S1_R1_001.bedGraph_pause_ratios_PROD.data",
                      col_names=c('tx_name', 'strand', 'c_1_301', 'coverage_c1_301'), delim=" ") %>%
    subset(select=-strand)
c_2_301 <- read_delim("Control_2_S2_R1_001.bedGraph_pause_ratios_PROD.data",
                      col_names=c('tx_name', 'strand', 'c_2_301', 'coverage_c2_301'), delim=" ") %>%
    subset(select=-strand)
c_3_301 <- read_delim("Control_3_S3_R1_001.bedGraph_pause_ratios_PROD.data",
                      col_names=c('tx_name', 'strand', 'c_3_301', 'coverage_c3_301'), delim=" ") %>%
    subset(select=-strand)
p_1_301 <- read_delim("Taf_1_S4_R1_001.bedGraph_pause_ratios_PROD.data",
                      col_names=c('tx_name', 'strand', 'p_1_301', 'coverage_p1_301'), delim=" ") %>%
    subset(select=-strand)
p_2_301 <- read_delim("Taf_2_S5_R1_001.bedGraph_pause_ratios_PROD.data",
                      col_names=c('tx_name', 'strand', 'p_2_301', 'coverage_p2_301'), delim=" ") %>%
    subset(select=-strand)
p_3_301 <- read_delim("Taf_3_S6_R1_001.bedGraph_pause_ratios_PROD.data",
                      col_names=c('tx_name', 'strand', 'p_3_301', 'coverage_p3_301'), delim=" ") %>%
    subset(select=-strand)

ddt_301 <- list(c_1_301, c_2_301, c_3_301, p_1_301, p_2_301, p_3_301) %>% reduce(left_join, by="tx_name")
ddt_301$c_pause_mean <- rowMeans(ddt_301[c('c_1_301', 'c_2_301', 'c_3_301')], na.rm=TRUE)
ddt_301$p_pause_mean <- rowMeans(ddt_301[c('p_1_301', 'p_2_301', 'p_3_301')], na.rm=TRUE)
ddt_301$c_coverage_mean <- rowMeans(ddt_301[c('coverage_c1_301', 'coverage_c2_301', 'coverage_c3_301')], na.rm=FALSE)
ddt_301$p_coverage_mean <- rowMeans(ddt_301[c('coverage_p1_301', 'coverage_p2_301', 'coverage_p3_301')], na.rm=FALSE)

ddt_301 <- subset(ddt_301, !duplicated(tx_name))
ddt_301 <- ddt_301 %>% mutate(pause_diff = log2(c_pause_mean) - log2(p_pause_mean))
ddt_301$coverage_mean <- rowMeans(ddt_301[c('c_coverage_mean', 'p_coverage_mean')], na.rm=TRUE)

plt <- function(frame, prefix) {
    ## 301 scatter
    ggplot(data = frame, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
        geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
        theme_tufte() +
        scale_fill_gsea() +
        labs(x = "Control",
             y= "Treatment",
             title= "Change in Pausing Index") +
        scale_x_log10() + scale_y_log10()
    ggsave(
        str_c(prefix, "_scatter.pdf"),
        plot = last_plot(), device = "pdf")

    ## Calculate Significance
    k <- ks.test(frame$p_pause_mean, frame$c_pause_mean)
    p <- k$p.value
    annotation <- paste0("p = ", p)
    print(annotation)
    write(annotation, file = str_c(prefix, "_pvalue.txt"))
    
    ## 301 ecdf
    ggplot(data = frame) +
        stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
        annotate("text", x = -Inf, y = Inf, vjust = 1, hjust = 0, label = annotation) +
        theme_tufte() +
        scale_color_npg() +
        theme(legend.title=element_blank()) +
        labs(x = "Pausing Index",
             y= "Cumulative Distribution",
             title= "Cumulative Distribution") +
        scale_x_log10()
    ggsave(
        str_c(prefix, "_ecdf.pdf"),
        plot = last_plot(), device = "pdf")

    ## MA-ish plot
    ggplot(data = frame, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
        geom_bin2d(bins=75) + geom_hline(yintercept=0) +
        scale_fill_gsea() +
        theme_tufte() +
        labs(x = "Normalized Coverage Level",
             y= "Change in Pause Index",
             title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
        scale_x_log10() +
        ## scale_y_log10()
        ggsave(
            str_c(prefix, "_MA.pdf"),
            plot = last_plot(), device = "pdf")
}

plt(ddt_301, "allgenes")

xfer <- read_delim("drosophila_refseq_to_common_id.txt", delim=" ",
                   col_names = c("tx_name", "name"))
paused_genes <- read_delim("genelist-paused.txt", delim="\t")
prox_genes <- read_delim("genelist-Prox.txt", delim="\t")
dist_genes <- read_delim("genelist-Dist.txt", delim="\t")
gaga_genes <- read_delim("genome_matched_genes_NVNVMGNRMR.data", delim="\t",
                         col_names = c("name", "seq"))
inr_genes <- read_delim("genome_matched_genes_BBCABW.data", delim="\t",
                        col_names = c("name", "seq"))
tata_genes <- read_delim("genome_matched_genes_WWWW.data", delim="\t",
                         col_names = c("name", "seq"))

paused_ddt <- left_join(ddt_301, xfer) %>%
    left_join(paused_genes) %>%
    na.omit(cols = c(genename))
prox_ddt <- left_join(ddt_301, xfer) %>%
    left_join(prox_genes) %>%
    na.omit(cols = c(genename))
dist_ddt <- left_join(ddt_301, xfer) %>%
    left_join(dist_genes) %>%
    na.omit(cols = c(genename))
gaga_ddt <- left_join(ddt_301, xfer) %>%
    left_join(gaga_genes) %>%
    na.omit(cols = c(genename))
inr_ddt <- left_join(ddt_301, xfer) %>%
    left_join(inr_genes) %>%
    na.omit(cols = c(genename))
tata_ddt <- left_join(ddt_301, xfer) %>%
    left_join(tata_genes) %>%
    na.omit(cols = c(genename))

plt(paused_ddt, "paused")
plt(prox_ddt, "proximal")
plt(dist_ddt, "distal")

plt(gaga_ddt, "gaga")
plt(inr_ddt, "inr")
plt(tata_ddt, "tata")
