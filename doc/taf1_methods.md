# Methods Used for Analysis of PRO-Seq Data

## Processing of Sequencing Data

The initial processing of all sequencing data was performed using the NascentFlow Pipeline, a data processing pipeline written in the Groovy programming language. The code for this pipeline can be found at https://github.com/Dowell-Lab/Nascent-Flow, with analysis for this experiment performed at commit 3fe1b7. Data were mapped to the hg38 reference genome for human cells, and to the dm6 reference genome for Drosophila cells.


## Isoform Resolution

For the remainder of the analysis, only the maximally expressed isoform of each gene was considered. The maximally expressed isoform was determined by calculating the RPKM normalized expression over each isoform and selecting the one with the maximum RPKM expression.

## Pause Index Analysis

Pause indices were calculated using a fixed-window approach. The region from -30 to +300 base pairs around the annotated transcription start site (TSS) was defined as the 'paused region', and the region from +301kb to the annotated polyA site was defined as the elongation region. Pause index was calculated as the ratio of length-normalized reads in the 'paused region' to length-normalized reads in the 'gene region'. Subsets of genes containing promoter elements were found by searching across the reference sequence of each gene, looking for  promoter elements in their expected positions relative to the TSS. The following motifs were used for each promoter element:
  - TATA-like: WWWW
  - Initiator: BBCABW
  - Motif Ten Element: CGANC....CGG
  - Downstream Promoter Element: RGWYVT
  - GAGA Element: NVNVMGNRMR (from Figure 3A in Tsai, 2016, Epigenetics/Chromatin)

## Metagene Analysis

Each gene in the isoform-resolved reference sequence was divided into a fixed number of bins, and the utility featurecounts (CITE) was used to determine the total counts in those regions. The mean count and standard deviation of that mean were calculated, and all bins were then plotted along with that standard deviation.

## Principal Component Analysis

Principal component analysis was performed using the standard `prcomp` function provided by the `sva` package for the R programming language. Batch effects from the different days replicates were generated on were corrected using the `removeBatchEffect` function provided by the `limma` package from the R programming language.

## Differential Expression Analysis

Differential expression analysis was performed using the `DESeq2` package for the R programming language. Counts were generated using the utility featurecounts. Initial analysis using counts across the full annotated gene showed significant skew indicating that the baseline assumptions of the differential expression model did not hold. To correct this, counts in the region from +500 of the TSS to -500 from the TES were used to obtain suitable model weights. Those model weights were then used when performing differential expression across the full gene, which corrected the observed skew.

## Gene Set Enrichment Analysis

Gene set enrichment analysis (GSEA) was performed with the Broad
Institute’s GSEA software on the GenePattern Server using the pre-ranked
module. Log₂ fold-change values were used as the rank metric for all
genes and compared against the Hallmark gene sets database for
enrichment.
