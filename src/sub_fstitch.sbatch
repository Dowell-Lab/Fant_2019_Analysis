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
BedDir=/scratch/Users/zama8258/processed_nascent/bedtools
OutDir=/scratch/Users/zama8258/processed_nascent/fstitch
bed="$BedDir"/"$file".trim.bedGraph
bed5="$BedDir"/"$file".trim.5.bedGraph
train="$OutDir"/"$file".train
out="$OutDir"/"$file".fstitch
trainfile=/scratch/Users/zama8258/pause_analysis_src/taf1knockdown_FSTrain_2.txt
FStitch=/scratch/Users/zama8258/fstitch-revisions/FastReadStitcher/src/FStitch

module load gcc/7.1.0
module load mpich/3.2.1

echo "Training FStitch Model (""$SECONDS"" sec)"
"$FStitch" train -cm 200 --strand + -np 32 -r "$bed" -t "$trainfile" -o "$train"
echo "Segmenting Pos Strand (""$SECONDS"" sec)"
"$FStitch" segment --strand + -np 32 -r "$bed" -w "$train" -o "$out".pos.annot.bed
echo "Segmenting Neg Strand (""$SECONDS"" sec)"
"$FStitch" segment --strand - -np 32 -r "$bed" -w "$train" -o "$out".neg.annot.bed

# echo "Combining Pos and Neg Strands (""$SECONDS"" sec)"
# module load bedtools/2.25.0

# cat "$out".pos.bedGraph \
		# 		"$out".neg.bedGraph \
		# 		> "$out".unsorted.bedGraph

# sortBed -i "$out".unsorted.bedGraph \
		# 				> "$OutDir"/annot/"$file".fstitch.annot.bedGraph

echo "Done (""$SECONDS"" sec)"
