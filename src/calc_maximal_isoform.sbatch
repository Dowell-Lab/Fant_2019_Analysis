#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=32gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zama8258@colorado.edu

## This script contains code to filter out the maximal isoform of all
## genes in a bedfile using bedops bedmap.

set -euxo pipefail

## Start timing the script
SECONDS=0

## This is a logrging function that will display the elapsed time of
## the script and how long the script has been executing for on the
## server. Useful for tracking execution time in various scripts.
function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

# Make sure we have the necessary modules on the cluster
if ! type -t bedtools
then module load bedtools
fi


##############################
## Variables we always need ##
##############################

# This is the file we're interested in analyzing
# InterestFile=$0 ## TODO	argument parsing
#InterestFile=/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/$srr.tri.BedGraph
RefSeq=/scratch/Shares/dowell/nascent/hg38/hg38_refseq.bed4

#####################################
## Generate	Temporary Files in RAM ##
#####################################

## Here we generate the temporary directory used to store the temp
## files we use during file processing. This is a fast process because
## we're only reading and writing to RAM, which gives us IO on the
## order of 60Gb/s
logr "Running ""$InterestFile"" in Production Mode"
TmpDir=$(mktemp -d)

PosRefFile="$TmpDir""/""$(uuidgen)"
NegRefFile="$TmpDir""/""$(uuidgen)"

RootName=$(basename "$InterestFile" .bam)
# PosStrandSums="$TmpDir""/""$(uuidgen)"
# NegStrandSums="$TmpDir""/""$(uuidgen)"

# Cache the sum files	because they take the longest to generate
PosStrandSums=/scratch/Users/zama8258/processed_nascent/scratch/"$RootName".pos.count
NegStrandSums=/scratch/Users/zama8258/processed_nascent/scratch/"$RootName".neg.count

PosCountsMapped=0
NegCountsMapped=0

PosScalingFactor=0
NegScalingFactor=0

PosStrandMillionsNormalizedSums="$TmpDir""/""$(uuidgen)"
NegStrandMillionsNormalizedSums="$TmpDir""/""$(uuidgen)"

PosStrandFpkmNormalizedSums="$TmpDir""/""$(uuidgen)"
NegStrandFpkmNormalizedSums="$TmpDir""/""$(uuidgen)"

PosStrandFinalSorted="$TmpDir""/""$(uuidgen)"
NegStrandFinalSorted="$TmpDir""/""$(uuidgen)"


# Clean up temp files on exit
function cleanup {
		rm -rf "$TmpDir"
		logr "Deleted temporary directory $TmpDir"
}
# Register the cleanup function to be called on the EXIT signal
trap cleanup EXIT

###############################
## Start Performing Analysis ##
###############################

## Here, we perform strand-specific filtering to account for `bedtools
## multicov` being unable to write out separate results for each
## strand. This strips down the RefSeq annotation to a bed4 format.
logr "Filtering RefSeq Annotations by Strand for ""$InterestFile"	&
awk -v OFS='\t' '{if ($4 > 0) print $1, $2, $3, $4}' "$RefSeq" \
		> "$PosRefFile" &
awk -v OFS='\t' '{if ($4 < 0) print $1, $2, $3, $4}' "$RefSeq" \
		> "$NegRefFile" &
wait

## Here, we make sure that our bam file is indexed and if not, perform
## the indexing
if [ ! -f "$RootName".bai ]; then
		logr "No Index File found. Generating One."
		module load samtools
		samtools index "$InterestFile"
fi

if [ ! -f "$PosStrandSums" ] || [ ! -f "$NegStrandSums" ]; then
		## Here, we construct a bedfile from our input BAM, in order to
		## normalize it by FPKM. We do this using bedtools multicov on the
		## sample in a strand specific manner. Note that bedtools multicov
		## will not account for duplicate reads by default, so if we want RPKM
		## instead of FPKM we need to add	the	`-D` flag
		logr "Calculating Strand-Specific Coverage for BAM" &
		bedtools multicov	-s -bams	"$InterestFile"	-bed "$PosRefFile" > "$PosStrandSums" &
		bedtools multicov	-s -bams	"$InterestFile"	-bed "$NegRefFile" > "$NegStrandSums"	&
		wait
fi


## Next, we need to find the sum of all reads mapped in our regions of
## interest, for FPKM. We do this by adding up every mapped count in
## the 5th column of our bedfile from bedtools multicov
read -r PosCountsMapped <<< "$(awk	-F '\t' -v OFS='\t'	'{sum += $5} END {print sum}' "$PosStrandSums")"
read -r NegCountsMapped <<< "$(awk	-F '\t' -v OFS='\t'	'{sum += $5} END {print sum}' "$NegStrandSums")"

## Then, we generate a scaling factor for RPKM by dividing our mapped
## counts by 1 million
read -r PosScalingFactor <<< "$((PosCountsMapped / 1000000))"
read -r NegScalingFactor <<< "$((NegCountsMapped / 1000000))"

## We proceed to divide every count in our bedfile by the
## strand-specific scaling factor
awk	-F '\t' -v OFS='\t' -v norm="$PosScalingFactor" \
		'{print $1, $2, $3, $4, $5 / norm}' "$PosStrandSums" > "$PosStrandMillionsNormalizedSums"
awk	-F '\t' -v OFS='\t' -v norm="$NegScalingFactor" \
		'{print $1, $2, $3, $4, $5 / norm}' "$NegStrandSums" > "$NegStrandMillionsNormalizedSums"

## We finally take the sums (now normalized for millions mapped), then
## divide those values by gene length.
awk -F '\t' -v OFS='\t'	'{print $1, $2, $3, $4, $5 / ($3 - $2)}' \
		"$PosStrandMillionsNormalizedSums" >	"$PosStrandFpkmNormalizedSums"
awk -F '\t' -v OFS='\t'	'{print $1, $2, $3, $4, $5 / ($3 - $2)}' \
		"$NegStrandMillionsNormalizedSums" >	"$NegStrandFpkmNormalizedSums"

## With the heavy lifting of normalization done, we can now easily
## take the maximum FPKM value for each isoform. Since isoforms have
## the same geneid, we can just sort in a descending order then unique
## by a specific column to get the maximal isoforms. After that, we
## sort back into	the normal bedfile format.
logr	"Taking Maximum Isoform and Resorting"
sort -rk5 -u "$PosStrandFpkmNormalizedSums" | sort -k 1,1 -k2,2n > "$PosStrandFinalSorted"
sort -rk5 -u "$NegStrandFpkmNormalizedSums" | sort -k 1,1 -k2,2n > "$NegStrandFinalSorted"
wait

## Here we perform a basic unit test to make sure that our output is
## valid. Simply, we count the lines in each output file and make sure
## that the number of lines is equivalent to the number of
logr "Verifying that Output is Valid"
PosLines=$(wc -l "$PosStrandFinalSorted")
NegLines=$(wc -l "$NegStrandFinalSorted")
read -r PosUniqLines <<< "$(awk '{print $4}' "$PosStrandFinalSorted" |
															 uniq -c | awk '{sum+=$1} END {print sum}')"
read -r NegUniqLines <<< "$(awk '{print $4}' "$NegStrandFinalSorted" |
															 uniq -c | awk '{sum+=$1} END {print sum}')"

if [ "$PosLines" -eq "$PosUniqLines" ] && [ "$NegLines" -eq "$NegUniqLines" ]; then
		logr "Passed tests:"
else
		logr "Failed tests:"
fi
logr "Pos Lines: ""$PosLines"", and Uniq Lines: ""$PosUniqLines"
logr "Neg Lines: ""$NegLines"", and Uniq Lines: ""$NegUniqLines"

## If we made it this far, we are done
logr "Done executing"
