#!/bin/bash

echo "Submitting Metagene Jobs for TAF1 Knockdown"
fileDir=/scratch/Users/zama8258/processed_nascent/metagene
topGenes="$fileDir"/topGenes.bed
middleGenes="$fileDir"/middleGenes.bed
bottomGenes="$fileDir"/bottomGenes.bed
topOut=/scratch/Users/zama8258/pause_output/metagene_top.png
middleOut=/scratch/Users/zama8258/pause_output/metagene_middle.png
bottomOut=/scratch/Users/zama8258/pause_output/metagene_bottom.png

sbatch sub_metagene.sh "$topGenes" "$topOut"
sbatch sub_metagene.sh "$middleGenes" "$middleOut"
sbatch sub_metagene.sh "$bottomGenes" "$bottomOut"
