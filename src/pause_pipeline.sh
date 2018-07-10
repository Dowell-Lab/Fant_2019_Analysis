#!/bin/bash
#SBATCH -p short
#SBATCH -N 1 
#SBATCH -c 32
#SBATCH --mem=128gb
#SBATCH --mail-user=zama8258@colorado.edu

CoreScript=./calc_pausing_indices_core2008.r
FixwinScript=./calc_pausing_indices_fixwin.sh
VizScript=./figures.r
SRR=SRR2084576

## Echo our scripts
echo $CoreScript
echo $FixwinScript
echo $SRR

## Actually run it all
module load bedtools
Rscript $CoreScript --srr=$SRR
bash $FixwinScript --pus=100 --pds=300 --gds=2000 --srr=$SRR
bash $FixwinScript --pus=100 --pds=300 --gds=5000 --srr=$SRR
Rscript $VizScript --srr=$SRR
