#!/bin/bash
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32gb
#SBATCH --mail-user=zama8258@colorado.edu

## TATA Box
bash /scratch/Users/zama8258/pause_analysis_src/filter_by_motif.sh \
		 --start=-50 --end=0 \
		 --fasta=/scratch/Users/zama8258/hg38.fa \
		 --motif=TATAAW \
		 --file=/scratch/Users/zama8258/NCBI_RefSeq_UCSC_RefSeq_hg38.bed
## Inr
bash /scratch/Users/zama8258/pause_analysis_src/filter_by_motif.sh \
		 --start=-10 --end=10 \
		 --fasta=/scratch/Users/zama8258/hg38.fa \
		 --motif=BBCABW \
		 --file=/scratch/Users/zama8258/NCBI_RefSeq_UCSC_RefSeq_hg38.bed
## MTE
bash /scratch/Users/zama8258/pause_analysis_src/filter_by_motif.sh \
		 --start=0 --end=100 \
		 --fasta=/scratch/Users/zama8258/hg38.fa \
		 --motif=CSARCSSAACGS \
		 --file=/scratch/Users/zama8258/NCBI_RefSeq_UCSC_RefSeq_hg38.bed
## DPE
bash /scratch/Users/zama8258/pause_analysis_src/filter_by_motif.sh \
		 --start=25 --end=40 \
		 --fasta=/scratch/Users/zama8258/hg38.fa \
		 --motif=RGWYVT \
		 --file=/scratch/Users/zama8258/NCBI_RefSeq_UCSC_RefSeq_hg38.bed

# Convert IUPAC Notation to GREP Compatible Regex TODO Implementation
# sed 's/W/[A|T]/g; s/R/[G|A]/g; s/Y/[T|C]/g; s/K/[G|T]/g; s/M/[A|C]/g; s/S/[G|C]/g; s/B/[G|T|C]/g; s/D/[G|A|T]/g; s/H/[A|C|T]/g; s/V/[G|C|A]/g; s/N/[A|G|C|T]/g'
