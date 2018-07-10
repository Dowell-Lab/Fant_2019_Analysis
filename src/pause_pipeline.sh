#!/bin/bash
#SBATCH -p short
#SBATCH -N 1 
#SBATCH -c 32
#SBATCH --mem=128gb
#SBATCH --mail-user=zama8258@colorado.edu

## Analysis Scripts
CoreScript=./calc_pausing_indices_core2008.r
FixwinScript=./calc_pausing_indices_fixwin.sh
## Visualization Scripts
VizScript=./figures.r
## SRR Number
SRR=SRR2084576

## Echo the things we're using
echo $CoreScript
echo $FixwinScript
echo $VizScript
echo $SRR

## Load Required Modules
module load bedtools
## Use CORE 2008 Method
Rscript $CoreScript --srr=$SRR
## Use CHEN 2015 Fixed Window Method with 2kbp Gene Body
bash $FixwinScript --pus=100 --pds=300 --gds=2000 --srr=$SRR
## Use CHEN 2015 Fixed Window Method with 5kbp Gene Body
bash $FixwinScript --pus=100 --pds=300 --gds=5000 --srr=$SRR
## Visualize Results
Rscript $VizScript --srr=$SRR
