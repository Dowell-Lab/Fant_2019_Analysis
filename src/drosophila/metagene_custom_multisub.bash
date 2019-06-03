#!/bin/bash

set -euo pipefail

echo "Submitting Metagene Jobs for TAF1 Knockdown"
fileDir=/scratch/Users/zama8258/taf1_drosophila_pro_seq/output/metagene
outDir="$fileDir"

topGenes="$fileDir"/topGenes.bed
middleGenes="$fileDir"/middleGenes.bed
bottomGenes="$fileDir"/bottomGenes.bed

topOut="$outDir"/metagene_top.pdf
middleOut="$outDir"/metagene_middle.pdf
bottomOut="$outDir"/metagene_bottom.pdf

quartile1Genes="$fileDir"/quartile1Genes.bed
quartile2Genes="$fileDir"/quartile2Genes.bed
quartile3Genes="$fileDir"/quartile3Genes.bed
quartile4Genes="$fileDir"/quartile4Genes.bed

quartile1Out="$outDir"/metagene_quartile1Genes.pdf
quartile2Out="$outDir"/metagene_quartile2Genes.pdf
quartile3Out="$outDir"/metagene_quartile3Genes.pdf
quartile4Out="$outDir"/metagene_quartile4Genes.pdf

sbatch metagene_custom.bash --regions="$topGenes" --outfile="$topOut"
sbatch metagene_custom.bash --regions="$middleGenes" --outfile="$middleOut"
sbatch metagene_custom.bash --regions="$bottomGenes" --outfile="$bottomOut"

sbatch metagene_custom.bash --regions="$quartile1Genes" --outfile="$quartile1Out"
sbatch metagene_custom.bash --regions="$quartile2Genes" --outfile="$quartile2Out"
sbatch metagene_custom.bash --regions="$quartile3Genes" --outfile="$quartile3Out"
sbatch metagene_custom.bash --regions="$quartile4Genes" --outfile="$quartile4Out"

echo "Submitted all Jobs Successfully"
