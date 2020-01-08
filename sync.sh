#!/bin/bash
# Sync to FIJI
echo "PUSHING"
# rclone copy -P --skip-links --bwlimit=10M \
# 	/home/zach/dowell_lab/pausing_meta_analysis/src/ \
# 	fiji:/scratch/Users/zama8258/pause_analysis_src
# echo "PULLING"
# rclone copy -P --dry-run --skip-links --bwlimit=10M \
# 	fiji:/scratch/Users/zama8258/pause_analysis_src \
# 	/home/zach/dowell_lab/pausing_meta_analysis/src/
rsync -Pzae "ssh" /home/zach/dowell_lab/pausing_meta_analysis/src/* \
	fiji2:/scratch/Users/zama8258/pause_analysis_src
rsync -Pzae "ssh" fiji2:/scratch/Users/zama8258/pause_analysis_src/* \
	/home/zach/dowell_lab/pausing_meta_analysis/src/
# rsync -Phe "ssh -T -o Compression=no -x" fiji:/scratch/Users/zama8258/pause_output/* out
