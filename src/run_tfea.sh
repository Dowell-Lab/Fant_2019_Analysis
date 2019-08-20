#!/bin/bash

set -euxo pipefail

outDir=/scratch/Users/zama8258/processed_nascent/tfea
bedDir=/scratch/Users/zama8258/processed_nascent_testing/mapped/bedgraphs
bamDir=/scratch/Users/zama8258/processed_nascent_testing/mapped/bams
controlBed1="$bedDir"/PO_1_S1_R1_001.bedGraph
controlBed2="$bedDir"/PO_2_S2_R1_001.bedGraph
treatmentBed1="$bedDir"/C413_1_S3_R1_001.bedGraph
treatmentBed2="$bedDir"/C413_2_S4_R1_001.bedGraph
controlBam1="$bamDir"/PO_1_S1_R1_001.sorted.bam
controlBam2="$bamDir"/PO_2_S2_R1_001.sorted.bam
treatmentBam1="$bamDir"/C413_1_S3_R1_001.sorted.bam
treatmentBam2="$bamDir"/C413_2_S4_R1_001.sorted.bam

TFEA --output "$outDir" \
		 --bed1 "$controlBed1" "$controlBed2" \
		 --bed2 "$treatmentBed1" "$treatmentBed2" \
		 --bam1 "$controlBam1" "$controlBam2" \
		 --bam2 "$treatmentBam1" "$treatmentBam2" \
		 --label1 control --label2 knockdown \
		 --genomefasta /scratch/Shares/dowell/genomes/hg38/hg38.fa \
		 --fimo_motifs /scratch/Users/zama8258/meme/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme \
		 --cpus 10 --mem 100gb \
		 --md --mdd --sbatch zama8258@colorado.edu
