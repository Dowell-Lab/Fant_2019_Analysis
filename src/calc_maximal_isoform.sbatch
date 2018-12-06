#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zama8258@colorado.edu

## This script contains code to filter out the maximal isoform of all
## genes in a bedfile using bedops bedmap.

## Start timing the script
SECONDS=0

## How many cores we use
NUM_CORES=8

## This is a logging function that will display the elapsed time of
## the script and how long the script has been executing for on the
## server. Useful for tracking execution time in various scripts.
function log {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

# Make sure we have the necessary modules on the cluster
if ! type -t bedtools
then module load bedtools
fi


# Variables we always need
DirPrefix=/scratch/Users/zama8258
InterestFile=$srr
#InterestFile=/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/$srr.tri.BedGraph
Infile=$DirPrefix/NCBI_RefSeq_UCSC_RefSeq_hg38.bed
RefSeq=/scratch/Shares/dowell/nascent/hg38/hg38_refseq.bed
OutFile=$DirPrefix/pause_output/"$srr"_max_isoforms.data

## Here we generate the temporary directory used to store the temp
## files we use during file processing. This is a fast process because
## we're only reading and writing to RAM, which gives us IO on the
## order of 60Gb/s
log "Running ""$srr"" in Production Mode"
TmpDir=$(mktemp -d)

PosRefFile="$TmpDir""/""$(uuidgen)"
NegRefFile="$TmpDir""/""$(uuidgen)"

PosStrandFile="$TmpDir""/""$(uuidgen)"
NegStrandFile="$TmpDir""/""$(uuidgen)"

PosStrandSums="$TmpDir""/""$(uuidgen)"
NegStrandSums="$TmpDir""/""$(uuidgen)"

PosStrandMeans="$TmpDir""/""$(uuidgen)"
NegStrandMeans="$TmpDir""/""$(uuidgen)"

PosStrandNormalizedSums="$TmpDir""/""$(uuidgen)"
NegStrandNormalizedSums="$TmpDir""/""$(uuidgen)"

PosStrandFinalSorted="$TmpDir""/""$(uuidgen)"
NegStrandFinalSorted="$TmpDir""/""$(uuidgen)"

# Clean up temp files on exit
function cleanup {
		rm -rf "$TmpDir"
		log "Deleted temporary directory $TmpDir"
}
# Register the cleanup function to be called on the EXIT signal
trap cleanup EXIT

## Here, we perform strand-specific filtering to account for `bedtools
## map` being unable to	write out separate results for each strand
log "Filtering RefSeq Annotations by Strand for ""$srr"
awk -v OFS='\t' -v pus="$pus" -v pds="$pds" \
		'{if ($6 == "+") print $1, $2, $2, $4, $5, $6}' "$RefSeq"\
		| sort -k1,1 -k2,2n > "$PosRefFile" &
awk -v OFS='\t' -v gus="$gus" -v gds="$gds" \
		'{if ($6 == "-") print $1, $2, $3, $4, $5, $6}' "$RefSeq" \
		| sort -k1,1 -k2,2n > "$NegRefFile" &
log "Separating bedfile by Strand"
awk -v OFS='\t' '{if ($4 > 0) print $1, $2, $3, $4}' "$InterestFile" \
		> "$PosStrandFile" &
awk -v OFS='\t' '{if ($4 < 0) print $1, $2, $3, $4}' "$InterestFile" \
		> "$NegStrandFile" &
wait

log "Calculating RefSeq Region Sums for ""$srr"
bedtools map -a "$PosRefFile" -b "$NegStrandFile" -c $(($NUM_CORES/2))  -o sum \
		| awk -F '\t' '($5 != "." && $5 != 0) ' > "$PosStrandSums" &
bedtools map -a "$NegRefFile" -b "$NegStrandFile" -c $(($NUM_CORES/2))  -o sum \
		| awk -F '\t' '($5 != "." && $5 != 0) ' > "$PosStrandSums" &

log "Done executing"
