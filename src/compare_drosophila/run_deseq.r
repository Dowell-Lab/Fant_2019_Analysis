## Script to run DESeq2
## Maintainer: Zachary Maas <zama8258@colorado.edu>

#######################
## Preliminary Setup ##
#######################

library("tidyverse")
library("DESeq2")
library("ggthemes")

## The directory our counts table is in.
setwd('/home/zach/dowell_lab/pausing_meta_analysis/out/compare_drosophila/')
# Number of replicates
num_reps <- 3

## Set up a function for plotting
## This is used later to do all the analysis + plotting
makefig <- function(deseqdata, fileprefix) {
    res <- DESeq2::results(deseqdata)
    table(res$padj < 0.05)
    ## Order by adjusted p-value
    res <- res[order(res$padj), ]
    ## Merge with normalized count data
    resdata <- merge(as.data.frame(res), as.data.frame(counts(deseqdata, normalized = TRUE)),
                     by = "row.names", sort = FALSE)
    names(resdata)[1] <- "Gene"
    head(resdata)
    ## Write results in the default format
    write.csv(resdata, file = paste0(fileprefix, "-diffexpr-results.csv"))

                                        # Generate the MA plot
    p_thresh <- 0.05
    p_str_signif <- str_c("p < ", p_thresh)
    p_str_not_signif <- str_c("p >= ", p_thresh)
    black <- "#455a64"
    red <- "#ff1744"
    res$padj <- res$padj %>% replace_na(1)
    res_df <- as_tibble(res) %>% mutate(signif = padj < p_thresh)
    ggplot(data = res_df) +
        geom_point(aes(x = baseMean, y = log2FoldChange, color = signif)) +
        scale_x_log10() + scale_color_manual(values = c(black, red),
                                             name="Significance",
                                             labels=c(p_str_not_signif,
                                                      p_str_signif)) +
        theme_tufte() +
        theme(legend.position=c(0.9,0.9)) +
        theme(legend.background = element_rect(colour = "#eceff1",
                                               fill = "white", linetype='solid')) +
        labs(title = "DESeq2 Differential Expression",
             x = "Log10 Mean of Normalized Counts",
             y = "Log2 Fold-Change") +
        ggsave(paste0(fileprefix, "-diffexpr-maplot.pdf"), width = 10, height = 5)

    ## Generate a tsv with results
    write.table(na.omit(as.data.frame(res[2])),
                file = paste0(fileprefix, "-genes_for_gsea.txt"), row.names = TRUE,
                quote = FALSE)
    write.table(res, file = paste0(fileprefix, "-DESeq.res.txt"), append = FALSE, sep = "\t")

    ## Generate preranked file for GSEA
    rnkdf <- tibble(gene = rownames(res), rnk = res$pvalue / sign(res$log2FoldChange)) %>%
        arrange(desc(rnk)) %>% drop_na()
    write.table(rnkdf, file = paste0(fileprefix, ".rnk"),
                append = FALSE, col.names = FALSE, row.names = FALSE,
                quote = FALSE, sep = "\t")

    return(res)
}

## To filter out only the genes we've used and convert to common ID
idConvert <- read_delim("drosophila_refseq_to_common_id.txt",
                        col_names = c("rowname", "common"), delim = " ")
## Experimental condition
condition <- factor(c(rep("treated", num_reps), rep("untreated", num_reps)),
                    levels=c("treated", "untreated"))

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
full_seqdata <- read.delim("counts_full.txt_without_header", stringsAsFactors = FALSE, header = TRUE,
                           row.names = 1)
## Then, we filter to only include the resolved isoforms
full_df <- as_tibble(rownames_to_column(full_seqdata)) %>% left_join(idConvert)

control_df <- full_df %>% subset(select = c("Chr", "Start", "End",
                                            "Strand", "Length", "common",
                                            "Control_1_S1_R1_001.sorted.bam",
                                            "Control_2_S2_R1_001.sorted.bam",
                                            "Control_3_S3_R1_001.sorted.bam",
                                            "LacI_plus_Cu1_S1_R1_001.sorted.bam",
                                            "LacI_plus_Cu2_S2_R1_001.sorted.bam",
                                            "LacI_plus_Cu3_S3_R1_001.sorted.bam"))
treatment_df <- full_df %>% subset(select = c("Chr", "Start", "End",
                                              "Strand", "Length", "common",
                                              "Taf_1_S4_R1_001.sorted.bam",
                                              "Taf_2_S5_R1_001.sorted.bam",
                                              "Taf_3_S6_R1_001.sorted.bam",
                                              "TAF1_plus_Cu2_S5_R1_001.sorted.bam",
                                              "TAF1_plus_Cu3_S6_R1_001.sorted.bam",
                                              "TAF1_plus_Cu_S4_R1_001.sorted.bam"))

control_dt <- data.frame(control_df)
rownames(control_dt) <- control_dt$common
control_dt$common <- NULL

## Now, actually perform DE-Seq
control_countdata <- control_dt[, 6:ncol(control_dt)]
control_countdata <- as.matrix(control_countdata)

control_coldata <- data.frame(row.names = colnames(control_countdata), condition)

dds <- DESeqDataSetFromMatrix(countData = control_countdata,
                              colData = control_coldata,
                              design = ~condition)
dds <- DESeq(dds)

ctrl <- makefig(dds, "control")

treatment_dt <- data.frame(treatment_df)
rownames(treatment_dt) <- treatment_dt$common
treatment_dt$common <- NULL

## Now, actually perform DE-Seq
treatment_countdata <- treatment_dt[, 6:ncol(treatment_dt)]
treatment_countdata <- as.matrix(treatment_countdata)

treatment_coldata <- data.frame(row.names = colnames(treatment_countdata), condition)

dds <- DESeqDataSetFromMatrix(countData = treatment_countdata,
                              colData = treatment_coldata,
                              design = ~condition)
dds <- DESeq(dds)

treat <- makefig(dds, "treatment")

