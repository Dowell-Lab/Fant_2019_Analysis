## Submit all cram files for processing to dreg bw

script=/scratch/Users/zama8258/pause_analysis_src/bedgraph_to_dreg_bw.bash

for cram in /scratch/Users/zama8258/taf1_pipeline_tfit/mapped/crams/*.cram; do
		[ -e "$cram" ] || continue
		echo "Submitting ""$cram"
		sbatch --export=cram_file="$cram" "$script"
done
