#!/bin/bash

echo "Submitting Multiple CalcIsoform Jobs"
dirPrefix=/scratch/Users/zama8258/taf1_drosophila_pro_seq/mapped/bams
c1="$dirPrefix"/Control_1_S1_R1_001.sorted.bam
c2="$dirPrefix"/Control_2_S2_R1_001.sorted.bam
c3="$dirPrefix"/Control_3_S3_R1_001.sorted.bam
p1="$dirPrefix"/Taf_1_S4_R1_001.sorted.bam
p2="$dirPrefix"/Taf_2_S5_R1_001.sorted.bam
p3="$dirPrefix"/Taf_3_S6_R1_001.sorted.bam
sbatch --export	InterestFile="$c1" --job-name="C1-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$c2" --job-name="C2-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$c3" --job-name="C3-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$p1" --job-name="P1-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$p2" --job-name="P2-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$p3" --job-name="P3-MaxIso" calc_maximal_isoform.sbatch
