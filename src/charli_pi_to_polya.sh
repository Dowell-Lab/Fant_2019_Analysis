#!/bin/bash
# Pausing Ratio Calculator
# Author: Zachary Maas <zama8258@colorado.edu>
# Licensed under GPLv3

# FIXME - This is currently our main development version of this
# script, and is under heavy development

#	TODO
# -	Combine all modes into 1 script using an argument
#	- Use	tempfiles instead	of a fixed directory
# - Add better error handling
#	-	Make it all faster

# Strict error checking
set -e
set -o nounset
set -o errexit

# Argument Parsing
# -pus/-pds = pause upstream/downstream
# -gds = gene upstream
# FIXME - add back in	-gus with new modes

usage()
{
		echo "calc_pausing_indices_fixwin.sh - a script for calculating fixed-width pausing indices"
		echo "Example:"
		echo "    ./calc_pausing_indices_fixwin.sh --pus=-100 --pds=300 --gds=2000 --srr=SRR2084556"
		echo "Usage:"
		echo "    -h/--help -- Display this help message."
		echo "    --pus     -- Pausing bases upstream"
		echo "    --pds     -- Pausing bases downstream"
		echo "    --gds     -- Gene bases downstream"
		echo "    --srr     -- SRR file to parse"
		exit 0
}

while [ "$1" != "" ]; do
    PARAM=$(echo "$1" | awk -F= '{print $1}')
		VALUE=$(echo "$1" | awk -F= '{print $2}')
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
gus=$(echo "$pds+1" | bc)

echo "Found parameters:"
echo "Up: $pus, Down: $pds, Gene Up: $gus, Gene Down: $gds, SRR: $srr"

# Make sure we have the necessary modules
if ! type -t bedtools 
then module load bedtools
fi

# Variables
DirPrefix=/scratch/Users/zama8258
TmpDir=charli_pi
Infile=$DirPrefix/NCBI_RefSeq_UCSC_RefSeq_hg19.bed
OutFile=$DirPrefix/pause_output/"$srr"_pause_ratios_$gds.data
OutGeneFile=$DirPrefix/$TmpDir/"$srr"_"$gds"_tss.bed
OutBodyFile=$DirPrefix/$TmpDir/"$srr"_"$gds"_body.bed
InterestFile=$srr
#InterestFile=/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/$srr.tri.BedGraph
InterestFilePos=$DirPrefix/$TmpDir/"$srr"_"$gds"_interest_pos.bed
InterestFileNeg=$DirPrefix/$TmpDir/"$srr"_"$gds"_interest_neg.bed

echo Prefiltering Reference Sequence...
awk -v OFS='\t' -v pus="$pus" -v pds="$pds" '{if ($6 == "+") print $1, $2-pus, $2+pds, $4, $5, $6; else print $1, $3-pus, $3+pds, $4, $5, $6}' $Infile \
		| sort -u -k1,1 -k2,2n | awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4}' > "$OutBodyFile" &
awk -v OFS='\t' -v gus="$gus" -v gds="$gds" '{if ($6 == "+") print $1, $2+gus, $3, $4, $5, $6; else print $1, $3-gus, $3, $4, $5, $6}' $Infile \
		| sort -u -k1,1 -k2,2n | awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4}' > "$OutGeneFile" &
wait

echo Splitting data file into pos and neg...
awk -v OFS='\t' '{if ($4 > 0) print $1, $2, $3, $4}' "$InterestFile" > "$InterestFilePos" &
awk -v OFS='\t' '{if ($4 < 0) print $1, $2, $3, $4}' "$InterestFile" > "$InterestFileNeg" &
wait

# Output Filenames
GeneOutPos=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_gene_pos.bed
GeneOutNeg=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_gene_neg.bed
BodyOutPos=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_body_pos.bed
BodyOutNeg=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_body_neg.bed

echo Finding Region Sums...
bedtools map -a "$OutGeneFile" -b "$InterestFilePos" -c 4 -o sum \
		| awk '($5 != "." && $5 != 0) ' > "$GeneOutPos" &
bedtools map -a "$OutBodyFile" -b "$InterestFilePos" -c 4 -o sum \
		| awk '($5 != "." && $5 != 0) ' > "$BodyOutPos" &
bedtools map -a "$OutGeneFile" -b "$InterestFileNeg" -c 4 -o sum \
		| awk '($5 != "." && $5 != 0) ' > "$GeneOutNeg" &
bedtools map -a "$OutBodyFile" -b "$InterestFileNeg" -c 4 -o sum \
		| awk '($5 != "." && $5 != 0) ' > "$BodyOutNeg" &
wait

CoverageOutPos=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_coverage_pos.bed
CoverageOutNeg=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_coverage_neg.bed

# FIXME - Finish Coverage Statistics
echo Calculating Coverage Statistics
awk 'NR==FNR{a[NR]=$5;next}{print $5+a[FNR]}' "$GeneOutPos" "$BodyOutPos" > "$CoverageOutPos"	&
awk 'NR==FNR{a[NR]=$5;next}{print $5+a[FNR]}' "$GeneOutNeg" "$BodyOutNeg" > "$CoverageOutNeg"	&
wait

FinalPos=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_final_pos.bed
FinalNeg=$DirPrefix/$TmpDir/"$srr"_"$gds"_out_final_neg.bed

# Add	Coverage to	Bodies
paste "$BodyOutPos" "$CoverageOutPos" > "$FinalPos"	&
paste "$BodyOutNeg" "$CoverageOutNeg" > "$FinalNeg" &
wait

echo Calculating Pausing Index
awk -F '\t' 'FNR==NR{a[$4]=$5; next} ($4 in a) {print $4,"+",$5/a[$4], $6}' "$GeneOutPos" "$FinalPos" > "$OutFile"
awk -F '\t' 'FNR==NR{a[$4]=$5; next} ($4 in a) {print $4,"-",$5/a[$4], $6}' "$GeneOutNeg" "$FinalNeg" >> "$OutFile"
echo "Done $srr $gds"
