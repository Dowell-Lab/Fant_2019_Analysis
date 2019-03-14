#!/bin/bash

set -euo pipefail

echo "Submitting Metagene Jobs for TAF1 Knockdown"
fileDir=/scratch/Users/zama8258/processed_nascent/metagene
topGenes="$fileDir"/topGenes.bed
middleGenes="$fileDir"/middleGenes.bed
bottomGenes="$fileDir"/bottomGenes.bed
topGenesAnti="$fileDir"/topGenesAnti.bed
middleGenesAnti="$fileDir"/middleGenesAnti.bed
bottomGenesAnti="$fileDir"/bottomGenesAnti.bed
topOut=/scratch/Users/zama8258/pause_output/metagene_top.png
middleOut=/scratch/Users/zama8258/pause_output/metagene_middle.png
bottomOut=/scratch/Users/zama8258/pause_output/metagene_bottom.png
topOutAnti=/scratch/Users/zama8258/pause_output/metagene_top_anti.png
middleOutAnti=/scratch/Users/zama8258/pause_output/metagene_middle_anti.png
bottomOutAnti=/scratch/Users/zama8258/pause_output/metagene_bottom_anti.png

sbatch sub_metagene.sh "$topGenes" "$topOut"
sbatch sub_metagene.sh "$middleGenes" "$middleOut"
sbatch sub_metagene.sh "$bottomGenes" "$bottomOut"
# sbatch sub_metagene.sh "$topGenesAnti" "$topOutAnti"
# sbatch sub_metagene.sh "$middleGenesAnti" "$middleOutAnti"
# sbatch sub_metagene.sh "$bottomGenesAnti" "$bottomOutAnti"
echo "Submitted all Jobs Successfully"
