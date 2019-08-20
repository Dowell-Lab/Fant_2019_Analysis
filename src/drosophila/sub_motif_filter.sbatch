#!/bin/bash
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --mail-user=zama8258@colorado.edu
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err

FASTA=/scratch/Shares/public/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa
REFSEQ=/scratch/Users/zama8258/dm6_refseq.bed

## TATA Box
bash /scratch/Users/zama8258/pause_analysis_src/drosophila/filter_by_motif.sh \
		 --start=-50 --end=25 \
		 --fasta="$FASTA" \
		 --motif=WWWW \
		 --file="$REFSEQ"
## Inr
bash /scratch/Users/zama8258/pause_analysis_src/drosophila/filter_by_motif.sh \
		 --start=-50 --end=50 \
		 --fasta="$FASTA" \
		 --motif=BBCABW \
		 --file="$REFSEQ"
## GAGA
bash /scratch/Users/zama8258/pause_analysis_src/drosophila/filter_by_motif.sh \
		 --start=-100 --end=50 \
		 --fasta="$FASTA" \
		 --motif=NVNVMGNRMR \
		 --file="$REFSEQ"
# Consensus: NVNVMGNRMR
# Taken from Tsai, 2016
## MTE
bash /scratch/Users/zama8258/pause_analysis_src/drosophila/filter_by_motif.sh \
		 --start=10 --end=100 \
		 --fasta="$FASTA" \
		 --motif=CGANC....CGG \
		 --file="$REFSEQ"
## DPE
bash /scratch/Users/zama8258/pause_analysis_src/drosophila/filter_by_motif.sh \
		 --start=10 --end=100 \
		 --fasta="$FASTA" \
		 --motif=RGWYVT \
		 --file="$REFSEQ"
## Pause Button
bash /scratch/Users/zama8258/pause_analysis_src/drosophila/filter_by_motif.sh \
		 --start=0 --end=100 \
		 --fasta="$FASTA" \
		 --motif=KCGRWCG \
		 --file="$REFSEQ"
## Pause Button is from Levine 2008

# Convert IUPAC Notation to GREP Compatible Regex TODO Implementation
# sed 's/W/[A|T]/g; s/R/[G|A]/g; s/Y/[T|C]/g; s/K/[G|T]/g; s/M/[A|C]/g;
# s/S/[G|C]/g; s/B/[G|T|C]/g; s/D/[G|A|T]/g; s/H/[A|C|T]/g; s/V/[G|C|A]/g; s/N/[A|G|C|T]/g'
