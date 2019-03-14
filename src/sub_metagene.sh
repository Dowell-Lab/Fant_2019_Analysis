#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=64gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zama8258@colorado.edu

SECONDS=0
script=/scratch/Users/zama8258/pause_analysis_src/run_metagene.r

echo "Starting..."

Rscript "$script" -i "$1" -o "$2"

echo "Done (""$SECONDS"" sec)"
