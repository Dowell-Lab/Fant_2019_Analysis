#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --mail-user=zama8258@colorado.edu

#	Assume only	utf-8
export LC_ALL=C

## Try -10kb
usage()
{
		echo "metagene_custom.bash - generate metagene plots for provided files"
		echo "Example:"
		echo "    ./metagene_custom --regions=regions.bed --outfile=metagene_regions.png"
		echo "Usage:"
		echo "    -h/--help -- Display this help message."
		echo "    --regions -- Bedgraph reference file to use"
		echo "    --outfile -- Output file, can be .png or .pdf"
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
				--regions)
						regionFile=$VALUE
						;;
				--outfile)
						imgOut=$VALUE
						;;
				*)
						echo "ERROR: unknown parameter \"$PARAM\""
						usage
						exit 1
						;;
		esac
		shift
done

function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

logr "Parsed Params: ""$regionFile"", ""$imgOut"

logr "Starting Analysis"
srcDir=/scratch/Users/zama8258/pause_analysis_src/metagene_fixed
bamDir=/scratch/Users/zama8258/processed_nascent_testing/mapped/bams
numRegions=100
# tmpdir=$(mktemp -d)
tmpdir=/scratch/Users/zama8258/processed_nascent/fixed_metagene/testing

NUM_CORES=8
safFile="$tmpdir"/"$(basename $imgOut .pdf)"_region_split.saf
countsSenseOut="$tmpdir"/metagene_counts_sense_"$(basename $imgOut .pdf)".txt
countsSenseFix="$tmpdir"/metagene_counts_sense_fix_"$(basename $imgOut .pdf)".txt
countsAntiSenseOut="$tmpdir"/metagene_counts_antisense_"$(basename $imgOut .pdf)".txt
countsAntiSenseFix="$tmpdir"/metagene_counts_antisense_fix_"$(basename $imgOut .pdf)".txt

if [ ! -f "$countsSenseFix" ] || [ ! -f "$countsAntiSenseFix" ]; then

		module load python/3.6.3
		logr "Generating Segmented SAF File"
		python3 "$srcDir"/bedgraph_split_for_metagene.py \
						-f "$regionFile" \
						-n "$numRegions" \
						-o "$safFile"

		logr "Changing Directories"
		module load subread
		cd "$bamDir" || exit
		logr "Building Sense Counts Table from SAF"
		featureCounts \
				-T "$NUM_CORES" \
				-s 1 \
				--fracOverlap 0.51 \
				-F 'SAF' \
				-a "$safFile" \
				-o "$countsSenseOut" \
				C413_1_S3_R1_001.sorted.bam \
				C413_2_S4_R1_001.sorted.bam \
				PO_1_S1_R1_001.sorted.bam \
				PO_2_S2_R1_001.sorted.bam

		logr "Building Antisense Counts Table from SAF"
		featureCounts \
				-T "$NUM_CORES" \
				-s 2 \
				--fracOverlap 0.51 \
				-F 'SAF' \
				-a "$safFile" \
				-o "$countsAntiSenseOut" \
				C413_1_S3_R1_001.sorted.bam \
				C413_2_S4_R1_001.sorted.bam \
				PO_1_S1_R1_001.sorted.bam \
				PO_2_S2_R1_001.sorted.bam

		logr "Fixing Counts File"
		tail -n+2 "$countsSenseOut" > "$countsSenseFix"
		tail -n+2 "$countsAntiSenseOut" > "$countsAntiSenseFix"

fi

logr "Generating Figure"
Rscript "$srcDir"/metagene_graph_custom.r \
				-s "$countsSenseFix" \
				-a "$countsAntiSenseFix" \
				-n "$numRegions" \
				-o "$imgOut"

logr "Done..."
