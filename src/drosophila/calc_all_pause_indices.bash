#!/bin/bash
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=64gb
#SBATCH --mail-user=zama8258@colorado.edu
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err

# Bed files
beddir=/scratch/Users/zama8258/taf1_drosophila_pro_seq/mapped/bedgraphs
cr1="$beddir"/Control_1_S1_R1_001.bedGraph
cr2="$beddir"/Control_2_S2_R1_001.bedGraph
cr3="$beddir"/Control_3_S3_R1_001.bedGraph
pr1="$beddir"/Taf_1_S4_R1_001.bedGraph
pr2="$beddir"/Taf_2_S5_R1_001.bedGraph
pr3="$beddir"/Taf_3_S6_R1_001.bedGraph

# FPKM files
fpkmdir=/scratch/Users/zama8258/taf1_drosophila_pro_seq/scratch/fpkm
cr1f="$fpkmdir"/Control_1_S1_R1_001.sorted.isoform_max.bed
cr2f="$fpkmdir"/Control_2_S2_R1_001.sorted.isoform_max.bed
cr3f="$fpkmdir"/Control_3_S3_R1_001.sorted.isoform_max.bed
pr1f="$fpkmdir"/Taf_1_S4_R1_001.sorted.isoform_max.bed
pr2f="$fpkmdir"/Taf_2_S5_R1_001.sorted.isoform_max.bed
pr3f="$fpkmdir"/Taf_3_S6_R1_001.sorted.isoform_max.bed

outdir=/scratch/Users/zama8258/taf1_drosophila_pro_seq/output/pausing
script=/scratch/Users/zama8258/pause_analysis_src/calculate_pause_index_to_polya.sh
upstream=-30
downstream=300
tag=PROD

bash "$script" \
		 --ref="$cr1f" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$cr1" &
bash "$script" \
		 --ref="$cr2f" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$cr2" &
bash "$script" \
		 --ref="$cr3f" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$cr3" &
bash "$script" \
		 --ref="$pr1f" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$pr1" &
bash "$script" \
		 --ref="$pr2f" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$pr2" &
bash "$script" \
		 --ref="$pr3f" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$pr3" &
wait

# ScriptDir=/scratch/Users/zama8258/pause_analysis_src
# python "$ScriptDir"/convert_isoform.py -l "$ScriptDir"/refseq_to_common_id.txt \
		# 			 -f "$StrandsMerged" \
		# 			 -o "$FPKMCommonID"
