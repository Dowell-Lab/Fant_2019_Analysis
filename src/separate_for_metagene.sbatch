#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --mail-user=zama8258@colorado.edu
# separate_for_metagene.bash --- Separates FPKM normalized calls for metagene
#
# Filename: separate_for_metagene.bash
# Description: Separates FPKM normalized calls for metagene generation
# Author: Student Zachary Maas <zama8258@colorado.edu>
# Maintainer: Student Zachary Maas <zama8258@colorado.edu>
# Created: Thu Feb 28 11:25:26 2019 (-0700)
#

# Commentary:
#
# This script should be fed in the fpkm files from the script
# calc_maximal_isoform.sbatch, and it will generate bedfiles to be fed
# into metagene for both the initiation region and the rest of the
# gene, after some preliminary filtering.
#

# Code:

set -euxo pipefail

# Data read-in and tmp files
srcBed=/scratch/Users/zama8258/processed_nascent/fpkm/C413_1_S3_R1_001.trim.sorted.isoform_max
scratch=/scratch/Users/zama8258/processed_nascent/metagene/
filtered="$scratch"/lenfilter.fpkm
sorted="$scratch"/sorted.fpkm
topGenes="$scratch"/topGenes.bed
bottomGenes="$scratch"/bottomGenes.bed
middleGenes="$scratch"/middleGenes.bed
lendiff="$scratch"/lendiff

# Variable setup
topAndBottomThreshold=500 # Take the top "n" genes off the top and bottom
minLengthThreshold=2000 # Minimum gene length that we require

# Filter off genes shorter than the threshold
awk -v OFS='\t' -v thresh="$minLengthThreshold" \
		'{if (sqrt($3-$2)^2 > thresh) print $0}' "$srcBed" > "$filtered"
wc -l "$srcBed" > "$lendiff"
wc -l "$filtered" >> "$lendiff"

# Take the top and bottom "threshold" genes as defined above
# First, sort in reverse order
sort -grk7,7 "$filtered" > "$sorted"

head -n"$topAndBottomThreshold" "$sorted" | \
		sort -k1,1 -k2,2 | cut -f1-6 > "$topGenes"
tail -n"$topAndBottomThreshold" "$sorted" | \
		sort -k1,1 -k2,2 | cut -f1-6 > "$bottomGenes"
head -n-"$topAndBottomThreshold" "$sorted" | \
		tail -n+"$topAndBottomThreshold" | cut -f1-6 > "$middleGenes"

#
# separate_for_metagene.bash ends here
