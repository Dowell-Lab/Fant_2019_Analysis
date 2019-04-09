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
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --mail-user=zama8258@colorado.edu

# Handle Errors Properly:
# -e - exits on the first error
# -u - fails when unset variables are called
# -o pipefail - fails when error in a pipe
set -euo pipefail

# Recently tweaked for better POSIX compliance.
function logr() {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

# Dynamically load bedtools if it isn't available
if ! type -t bedtools; then
		module load bedtools
fi

# Set the number of available cores
NUM_CORES=8

# TODO - Convert this to argument-parsing based for publication
# Variable setup
offset=500
negoffset=500
# TmpDir=/scratch/Users/zama8258/processed_nascent/scratch/features
TmpDir=$(mktemp -d)
BamDir=/scratch/Users/zama8258/processed_nascent_testing/mapped/bams
OutGenePosFile="$TmpDir""/gene_pos.bed"
OutGeneNegFile="$TmpDir""/gene_neg.bed"
GeneFileInr="$TmpDir""/all_genes_inr.bed"
GeneFileBody="$TmpDir""/all_genes_genebody.bed"
BaseDir=/scratch/Users/zama8258/processed_nascent
Infile="$BaseDir"/fpkm/C413_1_S3_R1_001.trim.sorted.isoform_max.bed

outFull="$BaseDir"/counts/counts_full.txt
outInr="$BaseDir"/counts/counts_initiaton.txt
outGenebody="$BaseDir"/counts/counts_genebody.txt

safFull="$BaseDir"/counts/full.saf
safInr="$BaseDir"/counts/initiation.saf
safGenebody="$BaseDir"/counts/genebody.saf

logr "Adjusting Coordinates"

logr "Performing Coordinate Adjustment about 5' End"
logr "NOTICE: \$negoffset will not be used in this mode."
awk -v OFS='\t' -v offset="$offset" -v negoffset="$negoffset" \
		'{if ($6 == "+") print $1, $2-offset, $2+offset, $4, $5, $6}' "$Infile" \
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
				> "$OutGenePosFile" &
awk -v OFS='\t' -v offset="$offset" -v negoffset="$negoffset" \
		'{if ($6 == "-") print $1, $3-offset, $3+offset, $4, $5, $6}' "$Infile" \
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
				> "$OutGeneNegFile" &
wait
## Merge Positive and Negative Strand Files and sort them.
logr "Merging Adjusted Coordinates"
cat	"$OutGenePosFile" "$OutGeneNegFile"	| sort -k1,1 -k2,2 > "$GeneFileInr"

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
## Merge Positive and Negative Strand Files and sort them.
logr "Merging Adjusted Coordinates"
cat	"$OutGenePosFile" "$OutGeneNegFile"	| sort -k1,1 -k2,2 > "$GeneFileBody"

## Generate the SAF Files
awk -v OFS='\t' '{print $4, $1, $2, $3, $6}' "$GeneFileInr" > "$safInr"
awk -v OFS='\t' '{print $4, $1, $2, $3, $6}' "$GeneFileBody" > "$safGenebody"
awk -v OFS='\t' '{print $4, $1, $2, $3, $6}' "$Infile" > "$safFull"

## Change directory because FeatureCounts is picky about running in
## the same directory as the bams.
logr "Changing Directories"
pushd "$BamDir" || exit
logr "Intersecting the Reference Sequence"

# Finally, run featurecounts to actually do the counting.
# Does this featurecounts need to count only 51% overlapping reads?

# Initation
logr "Performing Featurecounts (Inr)"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		-T "$NUM_CORES" \
		-s 1 \
		-F 'SAF' \
		-a "$safInr" \
		-o "$outInr" \
		C413_1_S3_R1_001.sorted.bam \
		C413_2_S4_R1_001.sorted.bam \
		PO_1_S1_R1_001.sorted.bam \
		PO_2_S2_R1_001.sorted.bam

tail -n +2 "$outInr" > "$outInr"_without_header

logr "Performing Featurecounts (Gene Body)"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		-T "$NUM_CORES" \
		-s 1 \
		-F 'SAF' \
		-a "$safGenebody" \
		-o "$outGenebody" \
		C413_1_S3_R1_001.sorted.bam \
		C413_2_S4_R1_001.sorted.bam \
		PO_1_S1_R1_001.sorted.bam \
		PO_2_S2_R1_001.sorted.bam

tail -n +2 "$outGenebody" > "$outGenebody"_without_header

logr "Performing Featurecounts (Full Gene)"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		-T "$NUM_CORES" \
		-s 1 \
		-F 'SAF' \
		-a "$safFull" \
		-o "$outFull" \
		C413_1_S3_R1_001.sorted.bam \
		C413_2_S4_R1_001.sorted.bam \
		PO_1_S1_R1_001.sorted.bam \
		PO_2_S2_R1_001.sorted.bam

tail -n +2 "$outFull" > "$outFull"_without_header

popd

logr "Done"

#
# gen_featurecounts.sbatch ends here
