#!/bin/bash
#SBATCH --job-name=tfit 
#SBATCH --mail-type=ALL			 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=zama8258@colorado.edu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=16:00:00			 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=500gb 			 # Memory limit
#SBATCH --output=/scratch/Users/zama8258/tfit_runs/taf1_2/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/tfit_runs/taf1_2/e_and_o/%x_%j.err

############ SET NUMBER OF THREADS FOR JOB TO USE PER NODE ###############

export OMP_NUM_THREADS=16

#################### LOAD MODULES #######################################

module load gcc/7.1.0
module load mpich/3.2.1

#################### SET VARIABLES ######################################

PROJECT=/scratch/Users/zama8258/tfit_runs/taf1_2
ROOTNAME=$1
TFIT=/scratch/Users/zama8258/tfit-revisions/src/Tfit
BEDPATH=/scratch/Users/zama8258/processed_nascent/bedtools
TSSFILE=/scratch/Shares/dowell/genomes/hg38/hg38_refseq.bed
# FSTITCH=/scratch/Shares/dowell/nascent/Allen2014/analysis/fstitch/combined/${ROOTNAME}.fstitch.bed
# TFITOUT=/scratch/Users/zama8258/tfit_runs/taf1
# CONFIGFILE=/scratch/Shares/dowell/Tfit/configs/config_orig.txt

#/scratch/Shares/dowell/nascent/Andrysik2017/tfit_train.bed
#/scratch/Shares/dowell/Tfit/hg38_TSS_UCSC_refseq_1000bprange.bed

####################### JOB INFO ########################################

printf "\nSample ID: $ROOTNAME"
printf "\nDirectory: $PROJECT"
printf "\nBedGraph: ${BEDPATH}/${ROOTNAME}.trim.rpkm.bedGraph"
printf "\nTfit Version: $TFIT"
printf "\nProject: $PROJECT"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

########################################################################

echo Running bidir module
date

$TFIT bidir \
			--bedgraph ${BEDPATH}/${ROOTNAME}.trim.rpkm.bedGraph \
			-N ${ROOTNAME} \
			--output ${PROJECT}/${ROOTNAME}.prelimbidir.bed \
			--logOut ${PROJECT}/logs/ \
			--bidirectionals $TSSFILE

# --segment ${FSTITCH}

# -k ${PROJECT}/FStitch_Beds_hg38/overcall/combined/${ROOTNAME}.FStitch.bed      Now an optional parameter

echo bidir module complete

### model module will attempt to make better mu calls on bidirectionals from prelim input and model RNAPII behavior
### mpi run will split the 

echo Running Model module
date

mpirun -np 4 $TFIT model \
			 --bedgraph ${BEDPATH}/${ROOTNAME}.trim.rpkm.bedGraph \
			 -N ${ROOTNAME} \
			 --output ${PROJECT}/model/${ROOTNAME}.bidir.train.bed \
			 --logOut ${PROJECT}/model/logs/ \
			 -s ${PROJECT}/${ROOTNAME}.prelimbidir.bed \
			 --bidirectionals $TSSFILE \
			 --threads 16

echo model module complete
echo Job completed
date
date
