#!env bash

echo "Submitting Multiple FStitch Jobs for TAF1 Knockdown"
c1=C413_1_S3_R1_001
c2=C413_2_S4_R1_001
p1=PO_1_S1_R1_001
p2=PO_2_S2_R1_001
sbatch --export	file="$c1" --job-name="C1-FS" sub_fstitch.sbatch
sbatch --export	file="$c2" --job-name="C2-FS" sub_fstitch.sbatch
sbatch --export	file="$p1" --job-name="P1-FS" sub_fstitch.sbatch
sbatch --export	file="$p2" --job-name="P2-FS" sub_fstitch.sbatch
