## Script to run DESeq2
## Maintainer: Zachary Maas <zama8258@colorado.edu>

#######################
## Preliminary Setup ##
#######################

library("tidyverse")
library("DESeq2")

## Set up a function for plotting
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
    ## Write results
    write.csv(resdata, file = paste0(fileprefix, "-diffexpr-results.csv"))

    ptsize <- 60

    ## png(paste0(fileprefix, "-diffexpr-hist.png"), 1500, 1000, pointsize = ptsize)
    ## hist(res$pvalue, breaks = 50, col = "grey")
    ## dev.off()

    ## png('MA_plot.png') plotMA(res) dev.off()

    png(paste0(fileprefix, "-diffexpr-maplot.png"), 3000, 2000, pointsize = ptsize)
    DESeq2::plotMA(res, main = "MA Plot")
    dev.off()

    ## Generate a tsv
    write.table(na.omit(as.data.frame(res[2])),
                file = paste0(fileprefix, "-genes_for_gsea.txt"), row.names = TRUE,
                quote = FALSE)
    write.table(res, file = paste0(fileprefix, "-DESeq.res.txt"), append = FALSE, sep = "\t")
}

## To filter out only the genes we've used and convert to common ID
idConvert <- read_delim("refseq_to_common_id.txt",
                        col_names = c("rowname", "common"), delim = " ")
## Experimental condition
condition <- factor(c(rep("treated", 2), rep("untreated", 2)), levels=c("untreated", "treated"))

#################################################
## Initial Processing for Correct Size Factors ##
#################################################

## Here, we use the genebody only and save appropriate plots

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
correct_seqdata <- read.delim("counts_genebody.txt_without_header", stringsAsFactors = FALSE, header = TRUE,
                              row.names = 1)
## Then, we filter to only include the resolved isoforms
correct_df <- as.tibble(rownames_to_column(correct_seqdata))
correct_dt <- data.frame(correct_df)
rownames(correct_dt) <- correct_dt$common
correct_dt$common <- NULL

## Now, actually perform DE-Seq
correct_countdata <- correct_dt[, 6:ncol(correct_dt)]
correct_countdata <- as.matrix(correct_countdata)
correct_countdata <- correct_countdata %>% subset(select = -c(Length))

correct_coldata <- data.frame(row.names = colnames(correct_countdata), condition)

cds <- DESeqDataSetFromMatrix(countData = correct_countdata, colData = correct_coldata, design = ~condition)
cds <- DESeq(cds)

## Make a figure for this data
makefig(cds, "genebody")

## Extract the size factors
sizes <- c(sizeFactors(cds))

###########################################################
## Then, use those size factors for the full set of data ##
###########################################################

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
full_seqdata <- read.delim("counts_fullgene.txt_without_header", stringsAsFactors = FALSE, header = TRUE,
                           row.names = 1)
## Then, we filter to only include the resolved isoforms
full_df <- as.tibble(rownames_to_column(full_seqdata))
## full_df <- inner_join(ids, full_df, by="rowname")
## full_df <- inner_join(full_df, idConvert, by = "rowname") %>% subset(select = -c(rowname)) %>%
## distinct(common, .keep_all = TRUE)
full_dt <- data.frame(full_df)
rownames(full_dt) <- full_dt$common
full_dt$common <- NULL

## Now, actually perform DE-Seq
full_countdata <- full_dt[, 6:ncol(full_dt)]
full_countdata <- as.matrix(full_countdata)
full_countdata <- full_countdata %>% subset(select = -c(Length))

full_coldata <- data.frame(row.names = colnames(full_countdata), condition)

dds <- DESeqDataSetFromMatrix(countData = full_countdata, colData = full_coldata, design = ~condition)
sizeFactors(dds) <- sizes
dds <- DESeq(dds)

makefig(dds, "fullgene")

######################################################
## Finally, use the initiation region only          ##
######################################################

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
initiation_seqdata <- read.delim("counts_initiation.txt_without_header", stringsAsFactors = FALSE, header = TRUE,
                                 row.names = 1)
## Then, we filter to only include the resolved isoforms
initiation_df <- as.tibble(rownames_to_column(initiation_seqdata))
## initiation_df <- inner_join(ids, initiation_df, by="rowname")
## initiation_df <- inner_join(initiation_df, idConvert, by = "rowname") %>% subset(select = -c(rowname)) %>%
## distinct(common, .keep_all = TRUE)
initiation_dt <- data.frame(initiation_df)
rownames(initiation_dt) <- initiation_dt$common
initiation_dt$common <- NULL

## Now, actually perform DE-Seq
initiation_countdata <- initiation_dt[, 6:ncol(initiation_dt)]
initiation_countdata <- as.matrix(initiation_countdata)
initiation_countdata <- initiation_countdata %>% subset(select = -c(Length))

initiation_coldata <- data.frame(row.names = colnames(initiation_countdata), condition)

eds <- DESeqDataSetFromMatrix(countData = initiation_countdata, colData = initiation_coldata, design = ~condition)
sizeFactors(eds) <- sizes
eds <- DESeq(eds)

makefig(eds, "initiation")
