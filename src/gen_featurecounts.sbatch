#!/bin/bash
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
# offset=-1000
# noffset=-1000
offset=0
noffset=0
TmpDir=/scratch/Users/zama8258/processed_nascent/scratch/features
BedHome=/scratch/Users/zama8258
OutGenePosFile="$TmpDir""/gene_pos.bed"
OutGeneNegFile="$TmpDir""/gene_neg.bed"
GeneFile="$TmpDir""/all_genes.bed"
Infile=/scratch/Users/zama8258/processed_nascent/fpkm/C413_1_S3_R1_001.trim.sorted.isoform_max.bed
out=counts_top500_fullgene.txt

logr "Adjusting Coordinates"
## Do coordinate adjustment using AWK, so we only count the	"gene body"
# awk -v OFS='\t' -v offset="$offset" -v noffset="$noffset" \
		# 		'{if ($6 == "+") print $1, $2+offset, $3-noffset, $4, $5, $6}' "$Infile" \
		# 		| sort -k1,1 -k2,2n | \
		# 		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
		# 				> "$OutGenePosFile" &
# awk -v OFS='\t' -v offset="$offset" -v noffset="$noffset" \
		# 		'{if ($6 == "-") print $1, $2+noffset, $3-offset, $4, $5, $6}' "$Infile" \
		# 		| sort -k1,1 -k2,2n | \
		# 		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
		# 				> "$OutGeneNegFile" &
# wait
awk -v OFS='\t' -v offset="$offset" -v noffset="$noffset" \
		'{if ($6 == "+") print $1, $2+offset, $2-offset, $4, $5, $6}' "$Infile" \
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
				> "$OutGenePosFile" &
awk -v OFS='\t' -v offset="$offset" -v noffset="$noffset" \
		'{if ($6 == "-") print $1, $2+offset, $3-offset, $4, $5, $6}' "$Infile" \
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4, $5, $6}' \
				> "$OutGeneNegFile" &
wait
logr "Merging Adjusted Coordinates"
cat	"$OutGenePosFile" "$OutGeneNegFile"	| sort -k1,1 -k2,2 > "$GeneFile"

logr "Changing Directories"
cd /scratch/Users/zama8258/processed_nascent_testing/mapped/bams || exit
logr "Intersecting the Reference Sequence"
## Filter out the set of genes we have selected as maximal isoforms.
## We can only use a single set of genes for this, so we must assume
## that all	samples have the same expressed set.
bedtools intersect -a "$BedHome"/hg38_reference.gtf \
				 -b	"$GeneFile" \
				 > "$BedHome"/hg38_reference_adjusted.gtf
logr "Performing Featurecounts"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts -T 32 -s 1 -a "$BedHome"/hg38_reference_adjusted.gtf -o "$out" C413_1_S3_R1_001.sorted.bam C413_2_S4_R1_001.sorted.bam PO_1_S1_R1_001.sorted.bam PO_2_S2_R1_001.sorted.bam

logr "Done"
