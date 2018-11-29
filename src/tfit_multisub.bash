#!env bash

echo "Submitting TFit Jobs for TAF1 Knockdown"
c1=C413_1_S3_R1_001
c2=C413_2_S4_R1_001
p1=PO_1_S1_R1_001
p2=PO_2_S2_R1_001
sbatch tfit_latest.sbatch "$c1"
sbatch tfit_latest.sbatch "$c2"
sbatch tfit_latest.sbatch "$p1"
sbatch tfit_latest.sbatch "$p2"
