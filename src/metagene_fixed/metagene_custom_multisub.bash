#!/bin/bash

set -euo pipefail

echo "Submitting Metagene Jobs for TAF1 Knockdown"
fileDir=/scratch/Users/zama8258/processed_nascent/fixed_metagene
outDir="$fileDir"

allGenes_tss="$fileDir"/allGenes_tss.bed
topGenes_tss="$fileDir"/topGenes_tss.bed
bottomGenes_tss="$fileDir"/bottomGenes_tss.bed
middleGenes_tss="$fileDir"/middleGenes_tss.bed
allGenes_body="$fileDir"/allGenes_body.bed
topGenes_body="$fileDir"/topGenes_body.bed
bottomGenes_body="$fileDir"/bottomGenes_body.bed
middleGenes_body="$fileDir"/middleGenes_body.bed
allGenes_tes="$fileDir"/allGenes_tes.bed
topGenes_tes="$fileDir"/topGenes_tes.bed
bottomGenes_tes="$fileDir"/bottomGenes_tes.bed
middleGenes_tes="$fileDir"/middleGenes_tes.bed

allOut_tss="$outDir"/metagene_fixed_tss_all.pdf
topOut_tss="$outDir"/metagene_fixed_tss_top.pdf
middleOut_tss="$outDir"/metagene_fixed_tss_middle.pdf
bottomOut_tss="$outDir"/metagene_fixed_tss_bottom.pdf
allOut_body="$outDir"/metagene_fixed_body_all.pdf
topOut_body="$outDir"/metagene_fixed_body_top.pdf
middleOut_body="$outDir"/metagene_fixed_body_middle.pdf
bottomOut_body="$outDir"/metagene_fixed_body_bottom.pdf
allOut_tes="$outDir"/metagene_fixed_tes_all.pdf
topOut_tes="$outDir"/metagene_fixed_tes_top.pdf
middleOut_tes="$outDir"/metagene_fixed_tes_middle.pdf
bottomOut_tes="$outDir"/metagene_fixed_tes_bottom.pdf

sbatch metagene_custom.bash --regions="$allGenes_tss" --outfile="$allOut_tss"
sbatch metagene_custom.bash --regions="$topGenes_tss" --outfile="$topOut_tss"
sbatch metagene_custom.bash --regions="$middleGenes_tss" --outfile="$middleOut_tss"
sbatch metagene_custom.bash --regions="$bottomGenes_tss" --outfile="$bottomOut_tss"
sbatch metagene_custom.bash --regions="$allGenes_body" --outfile="$allOut_body"
sbatch metagene_custom.bash --regions="$topGenes_body" --outfile="$topOut_body"
sbatch metagene_custom.bash --regions="$middleGenes_body" --outfile="$middleOut_body"
sbatch metagene_custom.bash --regions="$bottomGenes_body" --outfile="$bottomOut_body"
sbatch metagene_custom.bash --regions="$allGenes_tes" --outfile="$allOut_tes"
sbatch metagene_custom.bash --regions="$topGenes_tes" --outfile="$topOut_tes"
sbatch metagene_custom.bash --regions="$middleGenes_tes" --outfile="$middleOut_tes"
sbatch metagene_custom.bash --regions="$bottomGenes_tes" --outfile="$bottomOut_tes"

echo "Submitted all Jobs Successfully"
