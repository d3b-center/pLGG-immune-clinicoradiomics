#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Define directory and input files
data_dir="../../data"
count_file="${data_dir}/20230826_release-gene-counts-rsem-expected_count.collapsed.rds"
tpm_file="${data_dir}/20230826_release-gene-expression-rsem-tpm.collapsed.rds"
histology_file="${data_dir}/20230826_release.annotated_histologies.tsv"
tis_file="${data_dir}/tumor_inflammation_signatures.txt"

# script to run xcell on histology of interest
Rscript 00-run_xcell.R \
--mat $tpm_file \
--histology $histology_file \
--short_hist "LGAT" \
--broad_hist "Low-grade astrocytic tumor" \
--cg "low-grade glioma/astrocytoma" \
--output_dir "results/xcell_output"

# clustering of immune scores, combine xcell + clustering outputs and generate v-test scores
Rscript --vanilla 01-xcell_clustering.R \
--mat "results/xcell_output/xcell_scores.rds" \
--histology $histology_file \
--output_dir "results/clustering" \
--xcell_output_dir "results/xcell_output" \
--vtest_output_dir "results/vtest_analysis"

# generate heatmap of xcell scores annotated by molecular subtypes and clusters
Rscript --vanilla 02-xcell_heatmap.R \
--xcell_file "results/xcell_output/xcell_score_cluster.tsv" \
--plots_dir "plots/xcell_cluster_analysis"

# generate stats between xcell clusters vs clinical variables
Rscript --vanilla 03-xcell_cluster_analysis.R \
--xcell_file "results/xcell_output/xcell_score_cluster.tsv" \
--histology $histology_file \
--output_dir "results/xcell_cluster_analysis" \
--plots_dir "plots/xcell_cluster_analysis"

# Generate radarplots for v-test scores
Rscript --vanilla 04-vtest_plots.R \
--vtest_scores "results/vtest_analysis/vtest_scores_all.tsv" \
--plots_dir "plots/vtest_analysis"

# Compute TIS scores and generate plots
Rscript --vanilla 05-tis_analysis.R \
--count_file $count_file \
--xcell_file "results/xcell_output/xcell_score_cluster.tsv" \
--tis_genes $tis_file \
--output_dir "results/tis_analysis" \
--plots_dir "plots/tis_analysis"

# survival analysis using xCell-derived clusters
Rscript --vanilla 06-survival_analysis.R

# survival analysis using TIS high/low scores
Rscript --vanilla 06-survival_analysis_tis.R

# scatter plot of tumor purity and tumor fraction
Rscript --vanilla 07-tumor_purity_scatter.R \
--mat "results/xcell_output/xcell_scores.rds" \
--histology $histology_file \
--plots_dir "plots/tumor_purity_analysis" 

# differential pathway analysis on xCell-derived clusters
Rscript --vanilla 08-immuneSig_db_gsva.R \
--mat $tpm_file \
--xcell_file "results/xcell_output/xcell_score_cluster.tsv" \
--output_dir "results/gsva_analysis" \
--plots_dir "plots/gsva_analysis" 

# Compare imaging clusters with xCell-derived clusters
Rscript --vanilla 09-compare_classes.R
