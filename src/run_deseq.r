## Script to run DESeq2 Maintainter: Zachary Maas <zama8258@colorado.edu> Licensed
## under GPLv3

#######################
## Preliminary Setup ##
#######################

library("tidyverse")
library("DESeq2")

## To filter out only the genes we've used and convert to common ID
idConvert <- read_delim("/scratch/Users/zama8258/pause_analysis_src/refseq_to_common_id.txt",
                        col_names = c("rowname", "common"), delim = " ")
ids <- read_delim("C413_1_S3_R1_001.trim.bedGraph_pause_ratios_PROD.data",
                  col_names = c("rowname", "strand", "c_1_PROD", "coverage_c1_PROD"), delim = " ") %>%
    subset(select = rowname)

#################################################
## Initial Processing for Correct Size Factors ##
#################################################

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
correct_seqdata <- read.delim("counts_both_fix.txt", stringsAsFactors = FALSE, header = TRUE,
                              row.names = 1)
## Then, we filter to only include the resolved isoforms
correct_df <- as.tibble(rownames_to_column(correct_seqdata))
correct_df <- inner_join(ids, correct_df, by="rowname")
correct_df <- inner_join(correct_df, idConvert, by = "rowname") %>% subset(select = -c(rowname)) %>%
    distinct(common, .keep_all = TRUE)
correct_dt <- data.frame(correct_df)
rownames(correct_dt) <- correct_dt$common
correct_dt$common <- NULL

## Now, actually perform DE-Seq
correct_countdata <- correct_dt[, 6:ncol(correct_dt)]
correct_countdata <- as.matrix(correct_countdata)

condition <- factor(c(rep("treated", 2), rep("untreated", 2)), levels=c("untreated", "treated"))
correct_coldata <- data.frame(row.names = colnames(correct_countdata), condition)

cds <- DESeqDataSetFromMatrix(countData = correct_countdata, colData = correct_coldata, design = ~condition)
cds <- DESeq(cds)

## Extract the size factors
sizes <- c(sizeFactors(cds))

###########################################################
## Then, use those size factors for the full set of data ##
###########################################################

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
full_seqdata <- read.delim("counts_full_fix.txt", stringsAsFactors = FALSE, header = TRUE,
                           row.names = 1)
## Then, we filter to only include the resolved isoforms
full_df <- as.tibble(rownames_to_column(full_seqdata))
full_df <- inner_join(ids, full_df, by="rowname")
full_df <- inner_join(full_df, idConvert, by = "rowname") %>% subset(select = -c(rowname)) %>%
    distinct(common, .keep_all = TRUE)
full_dt <- data.frame(full_df)
rownames(full_dt) <- full_dt$common
full_dt$common <- NULL

## Now, actually perform DE-Seq
full_countdata <- full_dt[, 6:ncol(full_dt)]
full_countdata <- as.matrix(full_countdata)

condition <- factor(c(rep("treated", 2), rep("untreated", 2)), levels=c("untreated", "treated"))
full_coldata <- data.frame(row.names = colnames(full_countdata), condition)

dds <- DESeqDataSetFromMatrix(countData = full_countdata, colData = full_coldata, design = ~condition)
sizeFactors(dds) <- sizes
dds <- DESeq(dds)

######################################################
## With Proper Estimates, we can generate a MA Plot ##
######################################################

res <- results(dds)
table(res$padj < 0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), 
    by = "row.names", sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file = "diffexpr-results.csv")

png("diffexpr-hist.png", 1500, 1000, pointsize = 20)
hist(res$pvalue, breaks = 50, col = "grey")
dev.off()

## png('MA_plot.png') plotMA(res) dev.off()

png("diffexpr-maplot.png", 1500, 1000, pointsize = 20)
plotMA(res, main = "MA Plot")
dev.off()

## Volcano plot with 'significant' genes labeled
volcanoplot <- function(res, lfcthresh = 2, sigthresh = 0.05, main = "Volcano Plot", 
    legendpos = "bottomright", labelsig = TRUE, textcx = 1, ...) {
    with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = main, ...))
    with(subset(res, padj < sigthresh), points(log2FoldChange, -log10(pvalue), pch = 20, 
        col = "red", ...))
    with(subset(res, abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), 
        pch = 20, col = "orange", ...))
    with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), points(log2FoldChange, 
        -log10(pvalue), pch = 20, col = "green", ...))
    if (labelsig) {
        require(calibrate)
        with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), textxy(log2FoldChange, 
            -log10(pvalue), labs = Gene, cex = textcx, ...))
    }
    legend(legendpos, xjust = 1, yjust = 1, legend = c(paste("FDR<", sigthresh, sep = ""), 
        paste("|LogFC|>", lfcthresh, sep = ""), "both"), pch = 20, col = c("red", 
        "orange", "green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize = 20)
volcanoplot(resdata, lfcthresh = 1, sigthresh = 0.05, textcx = 0.8, xlim = c(-2.3, 
    2))
dev.off()

cat(rownames(subset(subset(res, pvalue < 0.01), log2FoldChange < 0)), sep = "\n")
cat(rownames(subset(subset(res, pvalue < 0.01), log2FoldChange > 0)), sep = "\n")

length(rownames(subset(subset(res, pvalue < 0.01), log2FoldChange < 0)))
length(rownames(subset(subset(res, pvalue < 0.01), log2FoldChange > 0)))
## Somme good genes: TOP: NM_001281749

## Generate a tsv
write.table(na.omit(as.data.frame(res[2])), file = "genes_for_gsea.txt", row.names = TRUE, 
    quote = FALSE)
write.table(res, file = "DESeq.res.txt", append = FALSE, sep = "\t")
