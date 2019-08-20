#!/bin/bash
# Sync to FIJI
rsync -Pzae "ssh" /home/zach/dowell_lab/pausing_meta_analysis/src/* \
	fiji2:/scratch/Users/zama8258/pause_analysis_src
rsync -Pzae "ssh" fiji:/scratch/Users/zama8258/pause_analysis_src/* \
	/home/zach/dowell_lab/pausing_meta_analysis/src/
# rsync -Phe "ssh -T -o Compression=no -x" fiji:/scratch/Users/zama8258/pause_output/* out
