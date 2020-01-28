# Pausing Meta Analysis

This project provides code used in the analysis of code from Fant 2019
(-- TITLE --). This analysis includes metagene analysis, pause index
analysis, and other miscellaneous analysis performed in the course of
developing this manuscript. Below, I have made a best effort to
describe the methods used in analysis and to group the files based on
the type of analysis that they belong to. If you have any questions
about specific code in this repository, please contact Zachary Maas at
<zama8258@colorado.edu>.

## Basic Project structure

- [src](./src) - scripts developed for data analysis
- [data](./data) - raw data used to perform analysis
- [doc](./doc) - documentation and preliminary reports

## Description of Individual Analyses

### Code Used In Publication
Code for pause index analysis:
 - <src/calculate_pause_index_to_polya.sh>
 - <src/gen_featurecounts.sbatch>
 - <src/gen_figures_final.r>
 - <src/gen_pause_index_figures.r>

Code for metagene analysis:
 - <src/bedgraph_split_for_metagene.py>
 - <src/metagene_custom.bash>
 - <src/metagene_graph_custom.r>
 - <src/normalizeBed.sh>

Code for differential expression analysis:
 - <src/run_deseq.r>
 - <src/run_deseq_compare_replicates.r>
 - <src/run_metagene.r>

Code for other analyses:
 - <src/calc_maximal_isoform.sbatch>
 - <src/run_pca_r.r>
 - <src/run_tfea.sh>
 - <src/separate_for_metagene.sbatch>
 - <src/make_moustache.r>
 - <src/bedgraph_to_dreg_bw.bash>
 - <src/hspa1b_fiveprime_read_ratio.r>
 - <src/hspa1b_select_reads_for_fiveprime_read_ratio.r>
 - <src/explore_lengths.r>
 - <src/check_long_gene_ends.r>
 - <src/convert_isoform.py>
 - <src/filter_by_motif.sh>
 - <src/filter_for_normalization.sbatch>

Supplemental files:
 - <src/refseq_to_common_id.txt>

#### Fixed window metagene code
 - <src/metagene_fixed/bedgraph_split_for_metagene.py>
 - <src/metagene_fixed/metagene_custom.bash>
 - <src/metagene_fixed/metagene_deeptools.bash>
 - <src/metagene_fixed/metagene_graph_custom.r>
 - <src/metagene_fixed/separate_for_metagene.sbatch>

#### Scripts for analysis of drosophila PRO-seq

Files here match their counterparts in the  - <src directory, adjusted for the drosophila data:>
 - <src/drosophila/calc_all_pause_indices.bash>
 - <src/drosophila/calc_maximal_isoform.sbatch>
 - <src/drosophila/filter_by_motif.sh>
 - <src/drosophila/gen_featurecounts.sbatch>
 - <src/drosophila/gen_pause_index_figures.r>
 - <src/drosophila/metagene_custom.bash>
 - <src/drosophila/metagene_graph_custom.r>
 - <src/drosophila/metaplot.bash>
 - <src/drosophila/metaplot_graph.r>
 - <src/drosophila/run_pca_r.r>
 - <src/drosophila/separate_for_metagene.sbatch>
 - <src/drosophila/sub_max_isoform.sbatch>
 - <src/drosophila/sub_motif_filter.sbatch>

Supplemental files:
 - <src/drosophila/genelist-Dist.bed>
 - <src/drosophila/genelist-Dist.txt>
 - <src/drosophila/genelist-Prox.bed>
 - <src/drosophila/genelist-Prox.txt>
 - <src/drosophila/genelist-paused.bed>
 - <src/drosophila_refseq_to_common_id.txt>

#### Scripts for analysis of metal shocked drosophila PRO-seq

Files here match their counterparts in the  - <src directory, adjusted for the drosophila metal shock data:>
 - <src/metal_drosophila/calc_all_pause_indices.bash>
 - <src/metal_drosophila/calc_maximal_isoform.sbatch>
 - <src/metal_drosophila/filter_by_motif.sh>
 - <src/metal_drosophila/gen_featurecounts.sbatch>
 - <src/metal_drosophila/gen_pause_index_figures.r>
 - <src/metal_drosophila/metagene_custom.bash>
 - <src/metal_drosophila/metagene_graph_custom.r>
 - <src/metal_drosophila/metaplot.bash>
 - <src/metal_drosophila/metaplot_graph.r>
 - <src/metal_drosophila/run_deseq.r>
 - <src/metal_drosophila/run_pca_r.r>
 - <src/metal_drosophila/separate_for_metagene.sbatch>
 - <src/metal_drosophila/sub_max_isoform.sbatch>
 - <src/metal_drosophila/sub_motif_filter.sbatch>

Supplemental files:
 - <src/metal_drosophila/genelist-Dist.bed>
 - <src/metal_drosophila/genelist-Dist.txt>
 - <src/metal_drosophila/genelist-Prox.bed>
 - <src/metal_drosophila/genelist-Prox.txt>
 - <src/metal_drosophila/genelist-paused.bed>

#### Scripts for comparison between drosophila experiments
 - <src/compare_drosophila/gen_featurecounts.sbatch>
 - <src/compare_drosophila/run_deseq.r>

### Scripts used to submit multiple analyses at once
 - <src/sub_fstitch.sbatch>
 - <src/sub_max_isoform.sbatch>
 - <src/sub_metagene.sh>
 - <src/sub_motif_filter.sbatch>
 - <src/sub_normalize_bed.sh>
 - <src/metagene_custom_multisub.bash>
 - <src/dreg_prep_multisub.sh>
 - <src/calc_all_pause_indices.bash>
 - <src/fstitch_multisub.bash>
 - <src/metagene_fixed/metagene_custom_multisub.bash>
 - <src/metagene_fixed/metagene_multisub.bash>
 - <src/drosophila/metagene_custom_multisub.bash>
 - <src/drosophila/metaplot_multisub.bash>
 - <src/metal_drosophila/metagene_custom_multisub.bash>
 - <src/metal_drosophila/metaplot_multisub.bash>

### Code not used in publication
 - <src/5prime_multisub.bash>
 - <src/separate_bedile.sh>
 - <src/setenv.bash>
 - <src/taf1knockdown_FSTrain.txt>
 - <src/taf1knockdown_FSTrain_2.txt>
 - <src/tagcounts.r>
 - <src/tfit-parser>
 - <src/tfit_latest.sbatch>
 - <src/tfit_multisub.bash>
 - <src/verify_isoform_differences.sh>
 - <src/pause_method_comparison.r>
 - <src/pause_pipeline.sh>
 - <src/prelim.r>
 - <src/charli_50bp_sliding_mean.sh>
 - <src/charli_figures.r>
 - <src/charli_figures_1kb.r>
 - <src/charli_pi_fixwin.sh>
 - <src/charli_pi_offset_tss.sh>
 - <src/charli_sub.sh>
 - <src/run_5primebed.sbatch>
 - <src/run_bedtools_test.sbatch>
 - <src/argparse_test.r>
 - <src/filter_counts_table.py>
 - <src/bootstrap.r>
 - <src/calc_pausing_indices.r>
 - <src/calc_pausing_indices_core2008.r>
 - <src/calc_pausing_indices_fixwin.sh>
 - <src/calc_pausing_indices_fixwin_norm.sh>
 - <src/calc_pausing_indices_tfit.R>
 - <src/get_std_chroms.bash>
 - <src/metagene_deeptools.bash>
 - <src/metagene_multisub.bash>
 - <src/run_pca.py>
 - <src/gen_comparison_figures.sbatch>

## License

Copyright (C) 2019 Dowell Lab

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
