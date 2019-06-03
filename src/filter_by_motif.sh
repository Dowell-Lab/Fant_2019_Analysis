#!/bin/bash
#	Utility to filter genome by	genes containing a known motif
# Author: Zachary Maas <zama8258@colorado.edu>
# Licensed under GPLv3

# Strict error checking
set -e
# set -o nounset
set -o errexit
#	Assume only	utf-8
export LC_ALL=C
# Echo start time

# Argument Parsing
# -start/-end = pause upstream/downstream
# -motif = gene upstream
# FIXME - add back in	-gus with new modes

usage()
{
		echo "filter_by_motif.sh - a script for filtering motifs"
		echo "Example:"
		echo "    ./filter_by_motif.sh --start=12 --end=30 --motif=TATAW --file=SRR0000000.bed"
		echo "Usage:"
		echo "    -h/--help -- Display this help message."
		echo "    --start   -- Bases upstream to search from"
		echo "    --end     -- Bases downstream to search to"
		echo "    --fasta   -- FASTA file to pull sequences from"
		echo "    --motif   -- Gene bases downstream"
		echo "    --file    -- Annotation file to parse"
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
				--start)
						start=$VALUE
						;;
				--end)
						end=$VALUE
						;;
				--fasta)
						fasta=$VALUE
						;;
				--motif)
						motif=$VALUE
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

echo "[PARAM] Start: $start, End: $end, Motif: $motif, File: $file"

# Make sure we have the necessary modules on the cluster
if ! type -t bedtools 
then module load bedtools
fi

################################################################################
################################################################################

# Variables we always need
DirPrefix=/scratch/Users/zama8258
#InterestFile=/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/$file.tri.BedGraph
# Infile=$DirPrefix/NCBI_RefSeq_UCSC_RefSeq_hg38.bed
InterestFile=$file
name=$(basename "$fasta" .fa)
OutFile=$DirPrefix/pause_output/"$name"_matched_genes_$motif.data

#	During debugging, we write out all output to disk so that we can
#	examine it and see what's going on with our script changes. This is
#	not necessary during production.
testing=false

if $testing; then

		# This DEBUG mode writes files in order to examine them to find
		# out where problems originated. This has a moderate performance
		# hit...
		echo "[LOG] Running ""$name"" in Debug Mode"
		TmpDir=charli_pi

		StrandGenes=$DirPrefix/$TmpDir/"$name"_"$motif"_genes.bed

		Sequences=$DirPrefix/$TmpDir/"$name"_"$motif"_seq.bed

		MatchedSequences=$DirPrefix/$TmpDir/"$name"_"$motif"_matched.bed

else

		# When not debugging, read and write to /tmp/ directly in memory
		# since it gives us a nice performance boost.
		echo "[LOG] Running ""$file"" in Production Mode."
		TmpDir=$(mktemp -d)

		StrandGenes="$TmpDir""/""$(uuidgen)"

		Sequences="$TmpDir""/""$(uuidgen)"

		MatchedSequences="$TmpDir""/""$(uuidgen)"

		# Clean up temp files on exit
		function cleanup {
				rm -rf "$TmpDir"
				echo "[LOG] Deleted temporary directory $TmpDir"
		}
		# Register the cleanup function to be called on the EXIT signal
		trap cleanup EXIT

fi

################################################################################
################################################################################

#	Find the regions we want to	search on both strands.
# FIXME - Generalize for non-human genomes...
echo "[LOG] Prefiltering ""$file"
grep -e "^chr[0-9\\|X\\|Y]*\\s" "$InterestFile" | \
		awk -v OFS='\t' -v start="$start" -v end="$end" -v motif="$motif" \
				'{if ($6 == "+") print $1, $2+start, $2+end, $4, $5, $6; else print $1, $3-end, $3-start, $4, $5, $6}' \
				> "$StrandGenes" &
wait

#	Generate the FASTA Sequences we're interested in, flipping negative
#	strands so we only need 1 query. This changes all negative
#	genes to have the "wrong" order, so we shouldn't use this file for
#	correlating with bedfiles
echo "[LOG] Extracting FASTA Sequences from ""$file"
bedtools getfasta -name -tab -s -bed "$StrandGenes"	\
				 -fi	"$fasta" -fo "$Sequences"	&
wait
# TODO - Check this...
# bedtools getfasta -bed "$NegativeStrandGenes"	-fi	"$fasta" -name -tab \
		# 		| perl -lane 'print $F[0] . " " . scalar reverse $F[1]'	\
		# 					 > "$NegativeSequences"	&

sub () {
		sed 's/W/[A\\|T]/g; s/R/[G\\|A]/g; s/Y/[T\\|C]/g; s/K/[G\\|T]/g; s/M/[A\\|C]/g; s/S/[G\\|C]/g; s/B/[G\\|T\\|C]/g; s/D/[G\\|A\\|T]/g; s/H/[A\\|C\\|T]/g; s/V/[G\\|C\\|A]/g; s/N/[A\\|G\\|C\\|T]/g'
}

motifpattern="$(echo "$motif" | sub)"

echo "[LOG] Filtering FASTA Sequences matching ""$motifpattern" &
grep -i	"$motifpattern" "$Sequences" > "$MatchedSequences" &
wait

cat	"$MatchedSequences" > "$OutFile"
