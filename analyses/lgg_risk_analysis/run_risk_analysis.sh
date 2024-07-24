#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# define directory and data files
data_dir="../../data"
tpm_file="${data_dir}/20230826_release-gene-expression-rsem-tpm.collapsed_subset.rds"
histology_file="${data_dir}/20230826_release.annotated_histologies.tsv"
risk_file="input/RiskScores_Grouping.xlsx"
gtf_file="${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz"

# script to run ssgsea
Rscript --vanilla 01-ssgsea.R \
--mat $tpm_file \
--histology_file $histology_file \
--risk_file $risk_file \
--gtf_file $gtf_file \
--output_dir "results"

# script to run elastic net regression
Rscript --vanilla 02-elastic_net_risk.R \
--histology_file $histology_file \
--risk_file $risk_file \
--output_dir "results" \
--plots_dir "plots"

# script to run germline analysis
Rscript --vanilla 03-germline_risk.R \
--histology_file $histology_file \
--risk_file $risk_file \
--germline_file "input/186_GOI_variant_summary.txt" \
--output_dir "results" \
--plots_dir "plots"
