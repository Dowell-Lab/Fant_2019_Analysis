#!/bin/bash
# setenv.bash --- Set Environment Variables for Data Analysis
#
# Filename: setenv.bash
# Description: Sets Environment Variables for TAF1 Data Analysis
# Author: Student Zachary Maas <zama8258@colorado.edu>
# Created: Thu Mar 14 15:55:59 2019 (-0600)
#

# Commentary:
#
# This file contains the code required to set environment variables
# for all analysis scripts used by the analysis of the TAF1 dataset.
# These scripts are all built to use this common set of environment
# variables. This file is sourced in all relevant scripts, and should
# set up everything to work correctly.
#

# Code:

# Set up proper error checking
set -euo pipefail

# First, we need to set up the data directory that's going to be used
# by the Nascentflow Nextflow pipeline. We'll use this directory to do
# the rest of the analysis, so it's crucial that this variable is set.
export BaseDir=/scratch/Users/zama8258/processed_nascent_testing

# Then, we set the directories for our various mapped datafiles.
export BamDir="$BaseDir"/mapped/bams
export BedDir="$BaseDir"/mapped/bedgraphs

# Check that all variables are correctly set. This blurb can be added
# to other scripts at the top to ensure that environment variables are
# properly set before attempting to call anything.
if [ ! -d "$BaseDir" ]; then
		echo "Base Directory ""$BaseDir"" not found."
elif [ ! -d "$BamDir" ]; then
		echo "Bam Directory ""$BamDir"" not found."
elif [ ! -d "$BedDir" ]; then
		echo "Bed Directory ""$BedDir"" not found."
fi

# Then, if they don't already exist, create our output files.
mkdir -p "$BaseDir"/scratch/{featurecounts,fpkm,e_and_o}
mkdir -p "$BaseDir"/output/{pca,pausing,metagene,diffexpr}

# TODO - Should I expand these to have their own variables?

# Finally, define some useful helper functions that we will use
# throughout the pipeline.

# A Logging Function for timing when output completes
function logr() {
    printf "[%s]: $*\n" "$(date -d@$SECONDS -u +%H:%M:%S)"
}
export -f logr

# Check if a module is loaded and load if not. This assumes that lmod
# is installed on the system you are using. If you don't want to do
# this, just change the function and use your locally installed
# versions of the required binaries instead.
function modcheck() {
		if ! type -t "$1"; then
				module load "$1"
		fi
}

modcheck bedtools

#
# setenv.bash ends here
