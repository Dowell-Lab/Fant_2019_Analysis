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
scratch=/scratch/Users/zama8258/processed_nascent/fixed_metagene
filtered_tss="$scratch"/lenfilter_tss.fpkm
filtered_body="$scratch"/lenfilter_body.fpkm
filtered_tes="$scratch"/lenfilter_tes.fpkm
sorted_tss="$scratch"/sorted_tss.fpkm
sorted_body="$scratch"/sorted_body.fpkm
sorted_tes="$scratch"/sorted_tes.fpkm

allGenes_tss="$scratch"/allGenes_tss.bed
topGenes_tss="$scratch"/topGenes_tss.bed
bottomGenes_tss="$scratch"/bottomGenes_tss.bed
middleGenes_tss="$scratch"/middleGenes_tss.bed
allGenes_body="$scratch"/allGenes_body.bed
topGenes_body="$scratch"/topGenes_body.bed
bottomGenes_body="$scratch"/bottomGenes_body.bed
middleGenes_body="$scratch"/middleGenes_body.bed
allGenes_tes="$scratch"/allGenes_tes.bed
topGenes_tes="$scratch"/topGenes_tes.bed
bottomGenes_tes="$scratch"/bottomGenes_tes.bed
middleGenes_tes="$scratch"/middleGenes_tes.bed

# lendiff="$scratch"/lendiff

# Variable setup
topAndBottomThreshold=500 # Take the top "n" genes off the top and bottom
# minLengthThreshold=2000 # Minimum gene length that we require

# Load bedtools
module load bedtools

# Filter About TSS
echo "Separating TSS"
awk -v OFS='\t' -f <(cat - <<- "EOF"
BEGIN {}
{
if (sqrt($3-$2)^2 > 2000 && $2 > 4001) {
  if ($6 == "+") {
    print $1, $2 - 2000, $2 + 2000, $4, $5, $6, $7, $8
  } else {
    print $1, $3 - 2000, $3 + 2000, $4, $5, $6, $7, $8
  }
}
}
END {}
EOF
										) "$srcBed" > "$filtered_tss"
# Take the top and bottom "threshold" genes as defined above
# First, sort in reverse order
sort -grk7,7 "$filtered_tss" > "$sorted_tss"

sort -k1,1 -k2,2 "$filtered_tss" | cut -f1-6 | bedtools sort > "$allGenes_tss" &
head -n"$topAndBottomThreshold" "$sorted_tss" | \
		sort -k1,1 -k2,2 | cut -f1-6 | bedtools sort > "$topGenes_tss" &
tail -n"$topAndBottomThreshold" "$sorted_tss" | \
		sort -k1,1 -k2,2 | cut -f1-6 | bedtools sort > "$bottomGenes_tss" &
head -n-"$topAndBottomThreshold" "$sorted_tss" | \
		tail -n+"$topAndBottomThreshold" | cut -f1-6 | \
		bedtools sort > "$middleGenes_tss" &
wait

echo "Separating Body"
awk -v OFS='\t' -f <(cat - <<- "EOF"
BEGIN {}
{
if (sqrt($3-$2)^2 > 4001) {
  if ($6 == "+") {
    print $1, $2 + 2000, $3 - 2000, $4, $5, $6, $7, $8
  } else {
    print $1, $2 + 2000, $3 - 2000, $4, $5, $6, $7, $8
  }
}
}
END {}
EOF
										) "$srcBed" > "$filtered_body"
# Take the top and bottom "threshold" genes as defined above
# First, sort in reverse order
sort -grk7,7 "$filtered_body" > "$sorted_body"

sort -k1,1 -k2,2 "$filtered_body" | cut -f1-6 | bedtools sort > "$allGenes_body" &
head -n"$topAndBottomThreshold" "$sorted_body" | \
		sort -k1,1 -k2,2 | cut -f1-6 | bedtools sort > "$topGenes_body" &
tail -n"$topAndBottomThreshold" "$sorted_body" | \
		sort -k1,1 -k2,2 | cut -f1-6 | bedtools sort > "$bottomGenes_body" &
head -n-"$topAndBottomThreshold" "$sorted_body" | \
		tail -n+"$topAndBottomThreshold" | cut -f1-6 | \
		bedtools sort > "$middleGenes_body" &
wait

echo "Separating TES"
awk -v OFS='\t' -f <(cat - <<- "EOF"
BEGIN {}
{
if (sqrt($3-$2)^2 > 2000 && $2 > 4001) {
  if ($6 == "+") {
    print $1, $3 - 2000, $3 + 2000, $4, $5, $6, $7, $8
  } else {
    print $1, $2 - 2000, $2 + 2000, $4, $5, $6, $7, $8
  }
}
}
END {}
EOF
										) "$srcBed" > "$filtered_tes"
# Take the top and bottom "threshold" genes as defined above
# First, sort in reverse order
sort -grk7,7 "$filtered_tes" > "$sorted_tes"

sort -k1,1 -k2,2 "$filtered_tes" | cut -f1-6 | bedtools sort > "$allGenes_tes" &
head -n"$topAndBottomThreshold" "$sorted_tes" | \
		sort -k1,1 -k2,2 | cut -f1-6 | bedtools sort > "$topGenes_tes" &
tail -n"$topAndBottomThreshold" "$sorted_tes" | \
		sort -k1,1 -k2,2 | cut -f1-6 | bedtools sort > "$bottomGenes_tes" &
head -n-"$topAndBottomThreshold" "$sorted_tes" | \
		tail -n+"$topAndBottomThreshold" | cut -f1-6 | \
		bedtools sort > "$middleGenes_tes" &
wait

echo "Done"

#
# separate_for_metagene.bash ends here
