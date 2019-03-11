#!/bin/bash
# gen_featurecounts.sbatch --- Generate Featurecounts
#
# Filename: gen_featurecounts.sbatch
# Description: Make Featurcounts with Coordinate Modification
# Author: Student Zachary Maas <zama8258@colorado.edu>
# Maintainer: Student Zachary Maas <zama8258@colorado.edu>
#

# Commentary:
#
# This file contains code for generating counts tables using
# featureCounts, with builtin support for coordinate modification for
# modifying regions in the bedfile.
#

# Code:

#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=32gb
#SBATCH --mail-user=zama8258@colorado.edu

set -euo pipefail

function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

module load bedtools

## Variable setup
offset=0
negoffset=0
TmpDir=/scratch/Users/zama8258/processed_nascent/scratch/features
BedHome=/scratch/Users/zama8258
OutGenePosFile="$TmpDir""/gene_pos.bed"
OutGeneNegFile="$TmpDir""/gene_neg.bed"
GeneFile="$TmpDir""/all_genes.bed"
Infile=/scratch/Users/zama8258/processed_nascent/fpkm/C413_1_S3_R1_001.trim.sorted.isoform_max.bed
out=counts_top500_fullgene.txt

logr "Adjusting Coordinates"

about_5prime=1 # TODO Change this to a script parameter
# Currently there are 2 modes that this script supports. Depending on
# how the about_5prime flag is set, coordinate adjustment will either
# be performed using a fixed-window about the TSS. If the flag is not
# set, coordinate adjustment will be performed using both values about
# the 5' and 3' ends.

# USAGE:
# about_5prime = 1
# Coordinates are TSS - offset to TSS + offset
# about_5prime = 0
# Coordinates are TSS + offset to TES - offset

if [ $about_5prime = 1 ]; then
		logr "Performing Coordinate Adjustment about 5' End"
		logr "NOTICE: \$negoffset will not be used in this mode."
		awk -v OFS='\t' -v offset="$offset" -v negoffset="$negoffset" \
				'{if ($6 == "+") print $1, $2+offset, $2-offset, $4, $5, $6}' "$Infile" \
				| sort -k1,1 -k2,2n | \
				awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
						> "$OutGenePosFile" &
		awk -v OFS='\t' -v offset="$offset" -v negoffset="$negoffset" \
				'{if ($6 == "-") print $1, $2+offset, $3-offset, $4, $5, $6}' "$Infile" \
				| sort -k1,1 -k2,2n | \
				awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
						> "$OutGeneNegFile" &
		wait
else
		logr "Performing Coordinate Adjustment about Full Gene"
		# Do coordinate adjustment using AWK, so we only count the	"gene body"
		awk -v OFS='\t' -v offset="$offset" -v negoffset="$negoffset" \
				'{if ($6 == "+") print $1, $2+offset, $3-negoffset, $4, $5, $6}' "$Infile" \
				| sort -k1,1 -k2,2n | \
				awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
						> "$OutGenePosFile" &
		awk -v OFS='\t' -v offset="$offset" -v negoffset="$negoffset" \
				'{if ($6 == "-") print $1, $2+negoffset, $3-offset, $4, $5, $6}' "$Infile" \
				| sort -k1,1 -k2,2n | \
				awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
						> "$OutGeneNegFile" &
		wait
fi

## Merge Positive and Negative Strand Files and sort them.
logr "Merging Adjusted Coordinates"
cat	"$OutGenePosFile" "$OutGeneNegFile"	| sort -k1,1 -k2,2 > "$GeneFile"

## Change directory because FeatureCounts is picky about running in
## the same directory as the bams.
logr "Changing Directories"
pushd /scratch/Users/zama8258/processed_nascent_testing/mapped/bams || exit
logr "Intersecting the Reference Sequence"

## Filter out the set of genes we have selected as maximal isoforms.
## We can only use a single set of genes for this, so we must assume
## that all	samples have the same expressed set.
bedtools intersect -a "$BedHome"/hg38_reference.gtf \
				 -b	"$GeneFile" \
				 > "$BedHome"/hg38_reference_adjusted.gtf

# Finally, run featurecounts to actually do the counting.
logr "Performing Featurecounts"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		-T 32 \
		-s 1 \
		-a "$BedHome"/hg38_reference_adjusted.gtf \
		-o "$out" \
		C413_1_S3_R1_001.sorted.bam \
		C413_2_S4_R1_001.sorted.bam \
		PO_1_S1_R1_001.sorted.bam \
		PO_2_S2_R1_001.sorted.bam
popd

# Make a copy of the final counts without header for ingesting into R
# scripts.
tail -n +2 "$out" "$out"_without_header

logr "Done"

#
# gen_featurecounts.sbatch ends here
