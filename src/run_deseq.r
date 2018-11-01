## Script to run DESeq2
## Maintainter: Zachary Maas <zama8258@colorado.edu>
## Licensed under GPLv3

library("tidyverse")
library("DESeq2")


seqdata <- read.delim("counts_fix.txt", stringsAsFactors = FALSE, header=TRUE,
                     row.names = 1)

countdata <- seqdata[ ,6:ncol(seqdata)]
countdata <- as.matrix(countdata)

condition <- factor(c(rep("untreated", 2), rep("treated", 2)))
coldata <- data.frame(row.names=colnames(countdata), condition)

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

dds <- DESeq(dds)

res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

png("diffexpr-hist.png", 1500, 1000, pointsize=20)
hist(res$pvalue, breaks=50, col="grey")
dev.off()

## png("MA_plot.png")
## plotMA(res)
## dev.off()

png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
plotMA(res, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
    with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
    with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
    if (labelsig) {
        require(calibrate)
        with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
    }
    legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

## High Expression: idx <- c(4651, 1562, 59, 154)
## "NR_003051" "NR_146151" "NR_137294" "NM_018667"
## RMRP, RNA45SN3, RNR1, SMPD3

## Low Expression: idx <- c(11, 13, 2, 2127)
## "NM_001270520" "NM_130900" "NM_001144966" "NM_003528"
## DAAM1, RAET1L, NEDD4L, HIST2H2BE

## TOP:    95  245  327  339  421  438  481  763  942 1282 1591 1721
## "NR_002586"    "NR_002450"    "NR_002918"    "NR_003940"    "NM_020643"
## "NM_052846"    "NM_001001710" "NM_080701"    "NM_053050"    "NR_136301"
## "NR_031729"    "NM_020862"
## SNORA63  SNORD68 SNORA48 SNORD80
## EMILIN3 FAM166A TREX2 MRPL53 LOC101929748
## MIR1908 LRFN1

## BOTTOM: 45   94   98  182  190  257  311  441  476  937  971 1365
## "NM_006795"    "NR_045566"    "NM_001199954" "NM_006389"    "NM_015299"
## "NM_001122636" "NR_038387"    "NR_024367"    "NM_139072"    "NM_003366"
## "NM_000700"    "NM_021009"
## EHD1 BBS4 ACTG1 HYOU1 KHNYN
## GALNT9 LINC00886 LINC00111 DNER UQCRC2
## ANXA1 UBC

cat(rownames(subset(subset(res, pvalue < 0.001), log2FoldChange < 0)), sep='\n')
cat(rownames(subset(subset(res, pvalue < 0.001), log2FoldChange > 0)), sep='\n')
## cat(rownames(subset(res, pvalue < 0.01)), sep='\n')

length(rownames(subset(subset(res, pvalue < 0.01), log2FoldChange < 0)))
length(rownames(subset(subset(res, pvalue < 0.01), log2FoldChange > 0)))
## Somme good genes:
## TOP: NM_001281749

## Generate a tsv
write.table(na.omit(as.data.frame(res[2])), file="genes_for_gsea.txt", row.names=TRUE, quote=FALSE)
write.table(res, file = "DESeq.res.txt", append = FALSE, sep= "	" )
