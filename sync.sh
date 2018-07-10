#!/bin/bash
# Sync to FIJI
rsync -Pzhe "ssh" src/* fiji:/scratch/Users/zama8258/pause_analysis_src
rsync -Pzhe "ssh" fiji:/scratch/Users/zama8258/pause_output/* out
