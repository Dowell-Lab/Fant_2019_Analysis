#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=100gb
#SBATCH --mail-user=zama8258@colorado.edu
# bedgraph_to_dreg_bw.bash --- Convert bedgraph to dreg suitable input
#
# Filename: bedgraph_to_dreg_bw.bash
# Description: Generate DREG input 5' begraph files
# Author: Zachary Maas <zama8258@colorado.edu>
# Maintainer: Zachary Maas <zama8258@colorado.edu>
# Created: Tue Jan 28 09:19:36 2020 (-0700)
#

# Commentary:
#
# Generate DREG input 5' begraph files
#

# Code:

set -euxo pipefail

module load bedtools

echo "Creating BigWigs suitable as inputs to dREG"

export CRAM_REFERENCE=/scratch/Shares/dowell/genomes/hg38/hg38.genome
chrom_sizes=/scratch/Shares/dowell/genomes/hg38/hg38.chrom.sizes
name=$(basename "$cram_file" .sorted.cram)
bedGraphToBigWig=/scratch/Users/zama8258/Nascent-Flow/bin/bedGraphToBigWig

pushd /scratch/Users/zama8258/taf1_pipeline_tfit/temp || exit
bamToBed -i "$cram_file" | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
		awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1, $3-1,$3,$4,$5,$6}' \
				> "$name".dreg.bed
sortBed -i "$name".dreg.bed > "$name".dreg.sort.bed

echo "positive strand processed to bedGraph"

bedtools genomecov \
         -bg \
         -i "$name".dreg.sort.bed \
         -g "$chrom_sizes" \
         -strand + \
         > "$name".pos.bedGraph
sortBed \
    -i "$name".pos.bedGraph \
    > "$name".pos.sort.bedGraph

bedtools genomecov \
         -bg \
         -i "$name".dreg.sort.bed \
         -g "$chrom_sizes" \
         -strand - \
    | awk 'BEGIN{FS=OFS="\t"} {$4=-$4}1' > "$name".neg.bedGraph
sortBed \
    -i "$name".neg.bedGraph \
    > "$name".neg.sort.bedGraph

echo "negative strand processed to bedGraph"

"$bedGraphToBigWig" "$name".pos.sort.bedGraph "$chrom_sizes" "$name".pos.bw
"$bedGraphToBigWig" "$name".neg.sort.bedGraph "$chrom_sizes" "$name".neg.bw

popd || exit

echo "bedGraph to bigwig done"
#
# bedgraph_to_dreg_bw.bash ends here
