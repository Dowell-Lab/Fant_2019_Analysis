### run_metagene.r --- Generate Metagene Plots
##
## Filename: run_metagene.r
## Description: Generate Metagene Plots for TAF1 Knockdown
## Author: Student Zachary Maas <zama8258@colorado.edu>
## Maintainer: Student Zachary Maas <zama8258@colorado.edu>
## Created: Thu Jan 31 11:37:49 2019 (-0700)
##
######################################################################
##
### Commentary:
##
## This file contains code to run metagene analysis on the 4
## replicates from our TAF1 knockdown experimet. The code uses the R
## metagene package.
##
### Code:

library("argparse")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", action="store", dest="regions",
                    help="The regions to generate a metagene plot for.")
parser$add_argument("-o", "--output", action="store", dest="outfile",
                    help="The output file name.")
args <- parser$parse_args()

regions <- c(args$regions)
outfile <- args$outfile
print(paste0("Regions: ", regions))
print(paste0("Out: ", outfile))

library("metagene")

print("Finishing Setup")
p1 <- "/scratch/Users/zama8258/processed_nascent_testing/mapped/bams/PO_1_S1_R1_001.sorted.bam"
p2 <- "/scratch/Users/zama8258/processed_nascent_testing/mapped/bams/PO_2_S2_R1_001.sorted.bam"
c1 <- "/scratch/Users/zama8258/processed_nascent_testing/mapped/bams/C413_1_S3_R1_001.sorted.bam"
c2 <- "/scratch/Users/zama8258/processed_nascent_testing/mapped/bams/C413_2_S4_R1_001.sorted.bam"

bam_files = c(p1, p2, c1, c2)
## regions <- c("/scratch/Users/zama8258/processed_nascent/fpkm/metagene_regions.bed")
## regions <- c("/scratch/Users/zama8258/processed_nascent/fpkm/top500.bed")
design <- data.frame(Samples = c(p1, p2, c1, c2),
                     control = c(1,1,0,0), treatment = c(0,0,1,1))

print("Generating Metagene")
mg <- metagene$new(regions = regions,
                   bam_files = bam_files,
                   force_seqlevels = TRUE)

mg$produce_table(design = design, flip_regions = TRUE)

print("Plotting Metagene")
## mg$plot(title = "Metagene Plot (Top 500 Genes)")
mg$plot(title = "Metagene Plot")

library("ggplot2")
## ggsave("/scratch/Users/zama8258/pause_output/metagene_top500.png", width=10, height=5)
ggsave(outfile, width=10, height=5)
## "/scratch/Users/zama8258/pause_output/metagene_all.png"

######################################################################
### run_metagene.r ends here
