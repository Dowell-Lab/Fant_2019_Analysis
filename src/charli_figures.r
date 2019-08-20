suppressMessages(library("tidyverse"))
library(ggthemes)
library(ggsci)

setwd('/scratch/Users/zama8258/pause_output/testing')
## Genesets we're interested in
wellexpr <- read_delim("expressed_genes_filtertable.txt", col_names=c("name", 'tx_name'), delim="\t")
inr <- read_delim("hg38_matched_genes_BBCABW.data", col_names=c("name", 'sequence'), delim="\t")
mte <- read_delim("hg38_matched_genes_CGANC....CGG.data", col_names=c("name", 'sequence'), delim="\t")
dpe <- read_delim("hg38_matched_genes_RGWYVT.data", col_names=c("name", 'sequence'), delim="\t")
tata <- read_delim("hg38_matched_genes_WWWW.data", col_names=c("name", 'sequence'), delim="\t")

## Filter and subset data
c_1_301 <- read_delim("C413_1_S3_R1_001.trim.bedGraph_pause_ratios_TEST.data",
                      col_names=c('tx_name', 'strand', 'c_1_301', 'coverage_c1_301'), delim=" ") %>%
    subset(select=-strand)
c_2_301 <- read_delim("C413_2_S4_R1_001.trim.bedGraph_pause_ratios_TEST.data",
                      col_names=c('tx_name', 'strand', 'c_2_301', 'coverage_c2_301'), delim=" ") %>%
    subset(select=-strand)
p_1_301 <- read_delim("PO_1_S1_R1_001.trim.bedGraph_pause_ratios_TEST.data",
                      col_names=c('tx_name', 'strand', 'p_1_301', 'coverage_p1_301'), delim=" ") %>%
    subset(select=-strand)
p_2_301 <- read_delim("PO_2_S2_R1_001.trim.bedGraph_pause_ratios_TEST.data",
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
    annotation <- paste0(prefix, ": p = ", p)
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

## Full gene sets
plt(ddt_301, "allgenes")
plt(ddt_301_wellexpressed, "wellexpressed")

## Single element subsets
plt(ddt_301_inr, "inr")
plt(ddt_301_mte, "mte")
plt(ddt_301_dpe, "dpe")
plt(ddt_301_tata, "tata")

plt(ddt_301_wellexpressed_inr, "wellexpressed_inr")
plt(ddt_301_wellexpressed_mte, "wellexpressed_mte")
plt(ddt_301_wellexpressed_dpe, "wellexpressed_dpe")
plt(ddt_301_wellexpressed_tata, "wellexpressed_tata")

## Multiple element subsets
plt(ddt_301_inr_and_mte, "inr_and_mte")
plt(ddt_301_inr_and_dpe, "inr_and_dpe")
plt(ddt_301_inr_and_tata, "inr_and_tata")
plt(ddt_301_tata_and_mte, "tata_and_mte")
plt(ddt_301_tata_and_dpe, "tata_and_dpe")

plt(ddt_301_wellexpressed_inr_and_mte, "wellexpressed_inr_and_mte")
plt(ddt_301_wellexpressed_inr_and_dpe, "wellexpressed_inr_and_dpe")
plt(ddt_301_wellexpressed_inr_and_tata, "wellexpressed_inr_and_tata")
plt(ddt_301_wellexpressed_tata_and_mte, "wellexpressed_tata_and_mte")
plt(ddt_301_wellexpressed_tata_and_dpe, "wellexpressed_tata_and_dpe")
