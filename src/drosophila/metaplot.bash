#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --mail-user=zama8258@colorado.edu

set -eo pipefail

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
		echo "    --regions -- Reference file to use"
		echo "    --outfile -- Output directory"
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
srcDir=/scratch/Users/zama8258/pause_analysis_src
numRegions=140
# tmpdir=$(mktemp -d)
tmpdir=/scratch/Users/zama8258/taf1_drosophila_pro_seq/output/metaplot_proxdist
prefix=$(basename "$regionFile" .bed)

NUM_CORES=8
safFile="$tmpdir"/"$prefix"_region_split.saf
countsSenseOut="$tmpdir"/"$prefix"_metagene_counts_sense.txt
countsSenseFix="$tmpdir"/"$prefix"_metagene_counts_sense_fix.txt
countsAntiSenseOut="$tmpdir"/"$prefix"_metagene_counts_antisense.txt
countsAntiSenseFix="$tmpdir"/"$prefix"_metagene_counts_antisense_fix.txt

module load python/3.6.3
logr "Generating Segmented SAF File"
python3 "$srcDir"/bedgraph_split_for_metagene.py \
				-f "$regionFile" \
				-n "$numRegions" \
				-o "$safFile"

logr "Changing Directories"
cd /scratch/Users/zama8258/taf1_drosophila_pro_seq/mapped/bams || exit
logr "Building Sense Counts Table from SAF"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		-T "$NUM_CORES" \
		-s 1 \
		-O \
		-F 'SAF' \
		-a "$safFile" \
		-o "$countsSenseOut" \
		Control_1_S1_R1_001.sorted.bam \
		Control_2_S2_R1_001.sorted.bam \
		Control_3_S3_R1_001.sorted.bam \
		Taf_1_S4_R1_001.sorted.bam \
		Taf_2_S5_R1_001.sorted.bam \
		Taf_3_S6_R1_001.sorted.bam

logr "Building Antisense Counts Table from SAF"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		-T "$NUM_CORES" \
		-s 2 \
		-O \
		-F 'SAF' \
		-a "$safFile" \
		-o "$countsAntiSenseOut" \
		Control_1_S1_R1_001.sorted.bam \
		Control_2_S2_R1_001.sorted.bam \
		Control_3_S3_R1_001.sorted.bam \
		Taf_1_S4_R1_001.sorted.bam \
		Taf_2_S5_R1_001.sorted.bam \
		Taf_3_S6_R1_001.sorted.bam

logr "Fixing Counts File"
tail -n+2 "$countsSenseOut" > "$countsSenseFix"
tail -n+2 "$countsAntiSenseOut" > "$countsAntiSenseFix"

logr "Generating Figure"
Rscript "$srcDir"/drosophila/metagene_graph_custom.r \
				-s "$countsSenseFix" \
				-a "$countsAntiSenseFix" \
				-n "$numRegions" \
				-o "$imgOut"

logr "Done..."
