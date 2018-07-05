#!/bin/bash
# Pausing Ratio Calculator
# Author: Zachary Maas <zama8258@colorado.edu>
# Licensed under GPLv3

# Strict error checking
set -e

# Make sure we have the necessary modules
if ! type -t bedtools 
then module load bedtools
fi

# Variables
DirPrefix=/scratch/Users/zama8258
TmpDir=genefilter_test
Infile=$DirPrefix/NCBI_RefSeq_UCSC_RefSeq_hg19.bed
OutFile=$DirPrefix/$TmpDir/pause_ratios.data
OutGeneFile=$DirPrefix/$TmpDir/tss.bed 
OutBodyFile=$DirPrefix/$TmpDir/body.bed 
InterestFile=/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/SRR2084576.tri.BedGraph
InterestFilePos=$DirPrefix/$TmpDir/interest_pos.bed 
InterestFileNeg=$DirPrefix/$TmpDir/interest_neg.bed 

echo Prefiltering Reference Sequence...
awk -v OFS='\t' '{if ($6 == "+") print $1, $2-100, $2+300, $4, $5, $6; else print $1, $3-100, $3+300, $4, $5, $6}' $Infile \
		| sort -k1,1 -k2,2n > $OutBodyFile &
awk -v OFS='\t' '{if ($6 == "+") print $1, $2+301, $2+2000, $4, $5, $6; else print $1, $3+301, $3+2000, $4, $5, $6}' $Infile \
		| sort -k1,1 -k2,2n > $OutGeneFile &
wait

echo Splitting data file into pos and neg...
awk -v OFS='\t' '{if ($4 > 0) print $1, $2, $3, $4}' $InterestFile > $InterestFilePos &
awk -v OFS='\t' '{if ($4 < 0) print $1, $2, $3, $4}' $InterestFile > $InterestFileNeg &
wait

# Output Filenames
GeneOutPos=$DirPrefix/$TmpDir/out_gene_pos.bed
GeneOutNeg=$DirPrefix/$TmpDir/out_gene_neg.bed
BodyOutPos=$DirPrefix/$TmpDir/out_body_pos.bed
BodyOutNeg=$DirPrefix/$TmpDir/out_body_neg.bed

echo Finding Region Sums...
bedtools map -a $OutGeneFile -b $InterestFilePos -c 4 -o sum \
		| awk '($7 != "." && $7 != 0) ' > $GeneOutPos &
bedtools map -a $OutGeneFile -b $InterestFileNeg -c 4 -o sum \
		| awk '($7 != "." && $7 != 0) ' > $GeneOutNeg &
bedtools map -a $OutBodyFile -b $InterestFilePos -c 4 -o sum \
		| awk '($7 != "." && $7 != 0) ' > $BodyOutPos &
bedtools map -a $OutBodyFile -b $InterestFileNeg -c 4 -o sum \
		| awk '($7 != "." && $7 != 0) ' > $BodyOutNeg &
wait

echo Calculating Pausing Index
# TODO - Change this to properly calculate pause ratios... Need to determine files.
awk -F '\t' 'FNR==NR{a[$4]=$7; next} ($4 in a) {print $4,"+",$7/a[$4]}' $GeneOutPos $BodyOutPos > $OutFile
awk -F '\t' 'FNR==NR{a[$4]=$7; next} ($4 in a) {print $4,"-",$7/a[$4]}' $GeneOutNeg $BodyOutNeg >> $OutFile

