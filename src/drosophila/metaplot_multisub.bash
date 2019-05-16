#!/bin/bash

set -euo pipefail

echo "Submitting Metagene Jobs for TAF1 Knockdown"
outDir=/scratch/Users/zama8258/taf1_drosophila_pro_seq/output/metaplot_proxdist
fileDir=/scratch/Users/zama8258/pause_analysis_src

proxGenes="$fileDir"/drosophila/genelist-Dist.bed
distGenes="$fileDir"/drosophila/genelist-Prox.bed
pausedGenes="$fileDir"/drosophila/genelist-paused.bed
topGenes=/scratch/Users/zama8258/taf1_drosophila_pro_seq/output/metagene/topGenes.bed

proxOut="$outDir"/metagene_prox.png
distOut="$outDir"/metagene_dist.png
pausedOut="$outDir"/metagene_paused.png
topOut="$outDir"/metagene_top.png

# sbatch metaplot.bash --regions="$proxGenes" --outfile="$proxOut"
# sbatch metaplot.bash --regions="$distGenes" --outfile="$distOut"
# sbatch metaplot.bash --regions="$pausedGenes" --outfile="$pausedOut"
sbatch metaplot.bash --regions="$topGenes" --outfile="$topOut"

echo "Submitted all Jobs Successfully"
