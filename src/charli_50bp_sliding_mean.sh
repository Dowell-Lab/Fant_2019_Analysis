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
		echo "aaa.sh - a script for calculating fixed-width window genome means"
		echo "Example:"
		echo "    ./charli_50bp_sliding_mean.sh --win=100 --slide=50 --file=x.bedGraph"
		echo "Usage:"
		echo "    -h/--help -- Display this help message."
		echo "    --win     -- Window size"
		echo "    --slide   -- How much to slide"
		echo "    --file    -- FILE file to parse"
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
				--win)
						win=$VALUE
						;;
				--slide)
						slide=$VALUE
						;;
				--file)
						file=$VALUE
						;;
				*)
						echo "ERROR: unknown parameter \"$PARAM\""
						usage
						exit 1
						;;
		esac
		shift
done

echo "Found parameters:"
echo "Win: $win, Slide: $slide, File: $file"

# Make sure we have the necessary modules
if ! type -t bedtools
then module load bedtools
fi

if ! type -t igvtools
then module load igvtools/2.3.75
fi

# Variables
# Genome=/scratch/Users/zama8258/processed_nascent/hg38.genome
DirPrefix=/scratch/Users/zama8258
TmpDir=charli_pi
OutFile=$DirPrefix/pause_output/$file.windows.bedGraph
Wins=$DirPrefix/$TmpDir/windows.bed
FileDir=/scratch/Users/zama8258/processed_nascent/bedtools/
InterestFile=$FileDir/$file
InterestFilePos=$DirPrefix/$TmpDir/$file\_50_interest_pos.bedGraph
InterestFileNeg=$DirPrefix/$TmpDir/$file\_50_interest_neg.bedGraph

bedtools makewindows -g $Genome -w $win -s $slide >	$Wins

echo Splitting data file into pos and neg...
awk -v OFS='\t' '{if ($4 > 0) print $1, $2, $3, $4}' $InterestFile > $InterestFilePos &
awk -v OFS='\t' '{if ($4 < 0) print $1, $2, $3, $4}' $InterestFile > $InterestFileNeg &
wait

echo "Mapping Positive to $OutFile"
bedtools map -a	$Wins -b $InterestFilePos -c	4	-o mean \
		| awk '($4 != "." && $4 != 0) '  \
					> $OutFile
echo "Mapping Negative to $OutFile"
bedtools map -a	$Wins -b $InterestFileNeg -c	4	-o mean \
		| awk '($4 != "." && $4 != 0) ' \
					>> $OutFile

sort -u -k1,1 -k2,2n > $OutFile.sort.bedGraph

echo "Generating TDF"
/opt/igvtools/2.3.75/igvtools toTDF $OutFile.sort.bedGraph \
															$OutFile.tdf \
		/scratch/Shares/dowell/genomes/hg38.chrom.sizes
