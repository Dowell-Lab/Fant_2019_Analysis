#!/bin/bash

set -euo pipefail

echo "Submitting Metagene Jobs for TAF1 Knockdown"
fileDir=/scratch/Users/zama8258/taf1_drosophila_pro_seq/output/metagene
outDir="$fileDir"

topGenes="$fileDir"/topGenes.bed
middleGenes="$fileDir"/middleGenes.bed
bottomGenes="$fileDir"/bottomGenes.bed

topOut="$outDir"/metagene_top.png
middleOut="$outDir"/metagene_middle.png
bottomOut="$outDir"/metagene_bottom.png

# shortGenes="$fileDir"/shortGenes.bed
# mediumGenes="$fileDir"/mediumGenes.bed
# mediumLongGenes="$fileDir"/mediumLongGenes.bed
# longGenes="$fileDir"/longGenes.bed

# shortOut="$outDir"/metagene_shortGenes.png
# mediumOut="$outDir"/metagene_mediumGenes.png
# mediumLongOut="$outDir"/metagene_mediumLongGenes.png
# longOut="$outDir"/metagene_longGenes.png

sbatch metagene_custom.bash --regions="$topGenes" --outfile="$topOut"
sbatch metagene_custom.bash --regions="$middleGenes" --outfile="$middleOut"
sbatch metagene_custom.bash --regions="$bottomGenes" --outfile="$bottomOut"

# sbatch metagene_custom.bash --regions="$shortGenes" --outfile="$shortOut"
# sbatch metagene_custom.bash --regions="$mediumGenes" --outfile="$mediumOut"
# sbatch metagene_custom.bash --regions="$mediumLongGenes" --outfile="$mediumLongOut"
# sbatch metagene_custom.bash --regions="$longGenes" --outfile="$longOut"

echo "Submitted all Jobs Successfully"
