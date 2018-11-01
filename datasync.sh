#!/bin/bash
# Sync to FIJI
# rsync -Phe "ssh -T -o Compression=no -x" src/* fiji:/scratch/Users/zama8258/pause_analysis_src
rsync -Praze "ssh" fiji:/scratch/Users/zama8258/pause_output/*.png out
rsync -Praze "ssh" fiji:/scratch/Users/zama8258/pause_output/hg38_matched_genes_* out
rsync -Praze "ssh" fiji:/scratch/Users/zama8258/pause_output/counts* out
