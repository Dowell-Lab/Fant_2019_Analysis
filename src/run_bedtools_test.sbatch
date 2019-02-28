#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=32gb
#SBATCH --mail-user=zama8258@colorado.edu

function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

module load bedtools

## Variable setup
## TmpDir=/scratch/Users/zama8258/processed_nascent/scratch/features
Infile=/scratch/Users/zama8258/processed_nascent/samtools/C413_1_S3_R1_001.trim.sorted.bam

logr "Generating BedGraph"
bedtools genomecov -ibam "$Infile" -bg -g /scratch/Shares/dowell/nascent/hg38/hg38.genome \
				 >	/scratch/Users/zama8258/processed_nascent/bedtools/tst.bed

logr "Done"
