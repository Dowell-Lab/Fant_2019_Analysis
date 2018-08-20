#!/bin/bash
# Pausing Ratio Calculator
# Author: Zachary Maas <zama8258@colorado.edu>
# Licensed under GPLv3

# Strict error checking
set -e

# Argument Parsing
# -pus/-pds = pause upstream/downstream
# -gds = gene upstream

usage()
{
		echo "calc_pausing_indices_fixwin_norm.sh - a script for calculating normalized fixed-width pausing indices"
		echo "Example:"
		echo "    ./calc_pausing_indices_fixwin_norm.sh --pus=-100 --pds=300 --gds=2000 --srr=SRR2084576"
		echo "Usage:"
		echo "    -h/--help -- Display this help message."
		echo "    --pus     -- Pausing bases upstream"
		echo "    --pds     -- Pausing bases downstream"
		echo "    --gds     -- Gene bases downstream"
		echo "    --srr     -- SRR file to parse"
		exit 0
}

while [ "$1" != "" ]; do
    PARAM=$(echo $1 | awk -F= '{print $1}')
		VALUE=$(echo $1 | awk -F= '{print $2}')
		case $PARAM in
				-h | --help)
						usage
						exit
						;;
				--pus)
						pus=$VALUE
						;;
				--pds)
						pds=$VALUE
						;;
				--gds)
						gds=$VALUE
						;;
				--srr)
						srr=$VALUE
						;;
				*)
						echo "ERROR: unknown parameter \"$PARAM\""
						usage
						exit 1
						;;
		esac
		shift
done

# Set Gene Downstream
gus=$pus+1

echo "Found parameters:"
echo "Up: $pus, Down: $pds, Gene Up: $gus, Gene Down: $gds, SRR: $srr"


# Make sure we have the necessary modules
if ! type -t bedtools 
then module load bedtools
fi

# Variables
DirPrefix=/scratch/Users/zama8258
TmpDir=genefilter_test
Infile=$DirPrefix/NCBI_RefSeq_UCSC_RefSeq_hg19.bed
OutFile=$DirPrefix/pause_output/$srr\_pause_ratios_$gds\_norm.data
OutGeneFile=$DirPrefix/$TmpDir/tss.bed 
OutBodyFile=$DirPrefix/$TmpDir/body.bed 
InterestFile=/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/$srr.tri.BedGraph
InterestFilePos=$DirPrefix/$TmpDir/interest_pos.bed 
InterestFileNeg=$DirPrefix/$TmpDir/interest_neg.bed 

echo Prefiltering Reference Sequence...
awk -v OFS='\t' -v pus="$pus" -v pds="$pds" '{if ($6 == "+") print $1, $2-pus, $2+pds, $4, $5, $6; else print $1, $3-pus, $3+pds, $4, $5, $6}' $Infile \
		| sort -k1,1 -k2,2n > $OutBodyFile &
awk -v OFS='\t' -v gus="$gus" -v gds="$gds" '{if ($6 == "+") print $1, $2+gus, $2+gds, $4, $5, $6; else print $1, $3+gus, $3+gds, $4, $5, $6}' $Infile \
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
# TODO create this
awk -F '\t' -v gus="$gus" -v gds="$gds" -v pus="$pus" -v pds="$pds" 'FNR==NR{a[$4]=$7; next} ($4 in a) {print $4,"+",($7/($pus+$pds)/(a[$4]/($gds-$gus))}' $GeneOutPos $BodyOutPos > $OutFile
awk -F '\t' -v gus="$gus" -v gds="$gds" -v pus="$pus" -v pds="$pds" 'FNR==NR{a[$4]=$7; next} ($4 in a) {print $4,"-",($7/($pus+$pds)/(a[$4]/($gds-$gus))}' $GeneOutNeg $BodyOutNeg >> $OutFile

