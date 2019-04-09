#!/bin/bash
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=64gb
#SBATCH --mail-user=zama8258@colorado.edu
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err

cr1=/scratch/Users/zama8258/processed_nascent/bedtools/C413_1_S3_R1_001.trim.rpkm.bedGraph
cr2=/scratch/Users/zama8258/processed_nascent/bedtools/C413_2_S4_R1_001.trim.rpkm.bedGraph
pr1=/scratch/Users/zama8258/processed_nascent/bedtools/PO_1_S1_R1_001.trim.rpkm.bedGraph
pr2=/scratch/Users/zama8258/processed_nascent/bedtools/PO_2_S2_R1_001.trim.rpkm.bedGraph

echo $cr1
echo $cr2
echo $pr1
echo $pr2

bash /scratch/Users/zama8258/pause_analysis_src/calculate_pause_index_to_polya.sh \
		 --ref=/scratch/Users/zama8258/processed_nascent/fpkm/C413_1_S3_R1_001.trim.sorted.isoform_max.bed \
		 --pus=-30 \
		 --pds=300 \
		 --gds=PROD \
		 --tag=C413_1_S3_R1_001.trim.bedGraph &
bash /scratch/Users/zama8258/pause_analysis_src/calculate_pause_index_to_polya.sh \
		 --ref=/scratch/Users/zama8258/processed_nascent/fpkm/C413_2_S4_R1_001.trim.sorted.isoform_max.bed \
		 --pus=-30 \
		 --pds=300 \
		 --gds=PROD \
		 --tag=C413_2_S4_R1_001.trim.bedGraph &
bash /scratch/Users/zama8258/pause_analysis_src/calculate_pause_index_to_polya.sh \
		 --ref=/scratch/Users/zama8258/processed_nascent/fpkm/PO_1_S1_R1_001.trim.sorted.isoform_max.bed \
		 --pus=-30 \
		 --pds=300 \
		 --gds=PROD \
		 --tag=PO_1_S1_R1_001.trim.bedGraph &
bash /scratch/Users/zama8258/pause_analysis_src/calculate_pause_index_to_polya.sh \
		 --ref=/scratch/Users/zama8258/processed_nascent/fpkm/PO_2_S2_R1_001.trim.sorted.isoform_max.bed \
		 --pus=-30 \
		 --pds=300 \
		 --gds=PROD \
		 --tag=PO_2_S2_R1_001.trim.bedGraph &
wait

# ScriptDir=/scratch/Users/zama8258/pause_analysis_src
# python "$ScriptDir"/convert_isoform.py -l "$ScriptDir"/refseq_to_common_id.txt \
# 			 -f "$StrandsMerged" \
# 			 -o "$FPKMCommonID"
