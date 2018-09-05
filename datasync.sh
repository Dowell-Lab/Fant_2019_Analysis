#!/bin/bash
# Sync to FIJI
# rsync -Phe "ssh -T -o Compression=no -x" src/* fiji:/scratch/Users/zama8258/pause_analysis_src
rsync -Phe "ssh -T -o Compression=no -x" fiji:/scratch/Users/zama8258/pause_output/*.png out
