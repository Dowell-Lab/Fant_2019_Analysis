#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/processed_nascent/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=32gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zama8258@colorado.edu

SECONDS=0
hg38=/scratch/Users/zama8258/processed_nascent/hg38.genome
BamDir=/scratch/Users/zama8258/processed_nascent/samtools
OutDir=/scratch/Users/zama8258/processed_nascent/bedtools
bam="$BamDir"/"$file"
out="$OutDir"/"$file"

module load bedtools/2.25.0

echo "Generating 5' BedGraphs (""$SECONDS"" sec)"
bedtools genomecov -bg -5	-strand + -ibam "$bam" -g "$hg38" >	"$out".pos.5.bedGraph
bedtools genomecov -bg -5	-strand - -ibam "$bam" -g "$hg38" >	"$out".neg.5.bedGraph

echo "Combining Pos and Neg Strands (""$SECONDS"" sec)"

cat "$out".pos.5.bedGraph \
		"$out".neg.5.bedGraph \
		> "$out".unsorted.5.bedGraph

sortBed -i "$out".unsorted.5.bedGraph \
				> "$out".sorted.5.bedGraph

echo "Done (""$SECONDS"" sec)"
