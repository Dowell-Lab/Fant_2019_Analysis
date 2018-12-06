# Methods Used for Analysis of PRO-Seq Data

## Processing of Sequencing Data

Sequencing data was initially in the form of raw fastq files. Prior to analysis, the data were analyzed using FastQC (CITE) to determine the quality of the data prior to trimming. Trimming of the data was performed using BBDuk (CITE), keeping reads with a minimum length of 25bp, a kmer length of 23, referenced against a standard collection of adapters for bbmap (CITE). FastQC was used to check the quality after trimming.

The trimmed fastq data was then mapped to the hg38 reference sequence using hisat2 (CITE), and these mapped data were sorted and converted to BAM format using samtools (CITE).

After mapping was completed, the mapped data was converted to the Bedgraph format using deeptools (CITE) for normalized bedgraph files and bedtools (CITE) for non-normalized and 5' bedgraph files.

## Pause Index Analysis

Initial analysis of the processed data focused on calculating the pause index based on NCBI RefSeq gene annotations. A fixed-window approach was used, taking the region from -50 to +200 base pairs around the annotated transcription start site (TSS) as the 'paused region' and taking the region from +1kb to the annotated polyA site as the 'gene region'. Pause index was calculated as the ratio of normalized reads in the 'paused region' to normalized reads in the 'gene region'.

Analysis of these results revealed that naive use of RefSeq annotations had a variety of issues. Most significantly, using these annotations meant that no resolution of actively coding gene isoforms could be determined, and that the annotated TSS frequently did not line up with the true TSS observed in gene browser traces.

These two issues (which appear to exist in any analysis naively using RefSeq annotations), meant that further analysis (like differential expression analysis) could not be performed without an effective method to determine the true TSS.

An alternative method for determining the true TSS was constructed using two existing algorithms -- FStitch (CITE) and TFit (CITE). The FStitch model was manually trained on the processed data to determine the regions in the genome most likely to be actively transcribed genes. Training was performed using a set of 66 annotated regions of the genome.
