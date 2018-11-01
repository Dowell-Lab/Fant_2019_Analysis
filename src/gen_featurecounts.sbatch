#!/bin/bash
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=32gb
#SBATCH --mail-user=zama8258@colorado.edu

# cd /scratch/Users/zama8258/processed_nascent/samtools/ && \
		# 		ls -l && \
		# 		/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		# 				-p \
		# 				-t exon \
		# 				-g gene_id  \
		# 				-T 32 \
		# 				-a ../../hg38_reference.gtf \
		# 				-o counts.txt C413_1_S3_R1_001.trim.sorted.bam C413_2_S3_R1_001.trim.sorted.bam PO_1_S1_R1_001.trim.sorted.bam PO_2_S2_R1_001.trim.sorted.bam

cd /scratch/Users/zama8258/processed_nascent/samtools/ || exit
echo "Filtering hg38"
grep start_codon ../../hg38_reference.gtf | \
		awk -v OFS='\t' '{if ($7 == "+") print $1, $2, $3, $4, $4 + 500, $6, $7, $8, $9" "$10" "$11" "$12; else print $1, $2, $3, $4, $4 + 500, $6, $7, $8, $9" "$10" "$11" "$12}' | \
		sed "s/start_codon/gene/" > ../../hg38_reference_adjusted.gtf
echo "Featurecounts"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts -T 32 -p -t gene -g gene_id  -a ../../hg38_reference_adjusted.gtf -o counts.txt PO_1_S1_R1_001.trim.sorted.bam PO_2_S2_R1_001.trim.sorted.bam C413_1_S3_R1_001.trim.sorted.bam C413_2_S4_R1_001.trim.sorted.bam
