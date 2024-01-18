# Perform consensus clustering on expected count of disease of interest 
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ConsensusClusterPlus)
  library(lspline)
  library(factoextra)
  library(fpc)
  library(Compind)
  library(COINr)
  library(IntNMF)
})

# parse arguments
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "input matrix (.rds)"),
  make_option(c("--histology"), type = "character",
              help = "histology file (.tsv)"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory for clustering analyses"),
  make_option(c("--xcell_output_dir"), type = "character",
              help = "output directory for xcell analyses"),
  make_option(c("--vtest_output_dir"), type = "character",
              help = "output directory for v-test analyses")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
mat <- opt$mat
histology <- opt$histology
output_dir <- opt$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = F)
xcell_output_dir <- opt$xcell_output_dir
vtest_output_dir <- opt$vtest_output_dir
dir.create(vtest_output_dir, recursive = TRUE, showWarnings = F)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# source functions from ClusTarIDseq
ClusTarIDseq_module <- file.path(root_dir, "analyses", "ClusTarIDseq")
source(file.path(ClusTarIDseq_module, "UQ.R"))
source(file.path(ClusTarIDseq_module, "run_ccp.R"))
source(file.path(ClusTarIDseq_module, "get_cdf_datapoints.R"))
source(file.path(ClusTarIDseq_module, "lspline_clustering.R"))
source(file.path(ClusTarIDseq_module, "dbscan_clustering.R"))
source(file.path(ClusTarIDseq_module, "nbmclust_clustering.R"))
source(file.path(ClusTarIDseq_module, "intnmf_clustering.R"))
source(file.path(ClusTarIDseq_module, "final_composite_score.R"))
source(file.path(ClusTarIDseq_module, "compute_vtest.R"))

# read histology
histology_df <- read_tsv(histology)

# # run lspline clustering and return the result with 
# lspline_clustering(expr_mat = mat, 
#                    hist_file = histology,
#                    filter_expr = FALSE, 
#                    protein_coding_only = FALSE,
#                    feature_selection = "variance",
#                    var_prop = 100, 
#                    transformation_type = "none", 
#                    max_k = 10, 
#                    coef_cutoff = 0.1, 
#                    min_cluster_size_prop = 0.1, max_cluster_size_prop = 0.5, 
#                    compute_all_equal = TRUE,
#                    output_dir = file.path(output_dir, "ccp_output"))
# 
# # run dbscan
# dbscan_clustering(expr_mat = mat, 
#                   filter_expr = FALSE,
#                   protein_coding_only = FALSE,
#                   feature_selection = "variance",
#                   var_prop = 100, 
#                   transformation_type = "none",
#                   minpts_val = NULL, 
#                   output_dir = file.path(output_dir, "dbscan_output"))
# 
# # run intNMF
# intnmf_clustering(expr_mat = mat, 
#                   filter_expr = FALSE, 
#                   protein_coding_only = FALSE, 
#                   feature_selection = "variance",
#                   var_prop = 100,
#                   transformation_type = "none",
#                   max_k = 10,
#                   output_dir = file.path(output_dir, "intnmf_output"))
# 
# # run NB.mclust is not applicable because it works with integers and 
# # it will convert all values to 0 in this case
# final_composite_score(expr_mat = mat, 
#                       hist_file = histology, 
#                       filter_expr = FALSE, 
#                       protein_coding_only = FALSE, 
#                       feature_selection = "variance",
#                       var_prop = 100,
#                       transformation_type = "none",
#                       ccp_output = file.path(output_dir, "ccp_output", "ccp_optimal_clusters.tsv"), 
#                       nbmclust_output = NULL, 
#                       dbscan_output = file.path(output_dir, "dbscan_output", "dbscan_optimal_clusters.tsv"),
#                       intnmf_output = file.path(output_dir, "intnmf_output", "intnmf_optimal_clusters.tsv"),
#                       output_dir = file.path(output_dir, "final_score"))

# read in results file 
clustering_method <- file.path(output_dir, "final_score", "final_clustering_output.tsv") %>%
  read_tsv() %>%
  filter(rank == 1) %>%
  pull(clustering_method)

if(clustering_method == "IntNMF"){
  final_output <- read_tsv(file.path(output_dir, "intnmf_output", "intnmf_optimal_clusters.tsv"))
} else if(clustering_method == "dbscan") {
  final_output <- read_tsv(file.path(output_dir, "dbscan_output", "dbscan_optimal_clusters.tsv"))
} else {
  final_output <- read_tsv(file.path(output_dir, "ccp_output", "ccp_optimal_clusters.tsv"))
}

# compute v test 
xcell_score <- readRDS(mat)
xcell_score <- xcell_score %>%
  rownames_to_column("cell_type") %>%
  gather(Kids_First_Biospecimen_ID, fraction, -c(cell_type)) %>%
  as.data.frame()
xcell_score <- xcell_score %>%
  dplyr::left_join(final_output, by = c("Kids_First_Biospecimen_ID" = "sample"))  %>%
  dplyr::group_by(cluster_assigned, cell_type) %>%
  dplyr::mutate(cluster_feature_mean_score = mean(fraction)) %>% # mean of gene count per cluster 
  ungroup() %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(feature_mean_score = mean(fraction),
                feature_variance = var(fraction)) 

# combine with histology molecular subtype
xcell_score <- histology_df %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort_participant_id, molecular_subtype) %>%
  inner_join(xcell_score) 

# reduce molecular subtypes by collapsing to larger groups
xcell_score$molecular_subtype <- as.character(xcell_score$molecular_subtype)
xcell_score$molecular_subtype_summarized <- xcell_score$molecular_subtype
# Group all samples containing 'BRAF V600E' in the subtype string into an 'LGG, BRAF V600E' group
xcell_score$molecular_subtype_summarized[grep("BRAF V600E", xcell_score$molecular_subtype)] <- "LGG, BRAF V600E"
# Group all samples containing 'BRAF-KIAA1549' in the subtype string into an 'LGG, BRAF-KIAA1549' group
xcell_score$molecular_subtype_summarized[grep("BRAF-KIAA1549", xcell_score$molecular_subtype)] <- "LGG, BRAF-KIAA1549"
# Group all remaining samples containing 'NF1' in the subtype string into an 'LGG_NF1' group.
xcell_score$molecular_subtype_summarized[grep("NF1", xcell_score$molecular_subtype)] <- "LGG_NF1"
# Group all remaining samples containing 'CDKN2A/B' in the subtype string into an 'LGG_CDKN2A/B' group.
xcell_score$molecular_subtype_summarized[grep("CDKN2A/B", xcell_score$molecular_subtype)] <- "LGG_CDKN2A/B"
# Group all remaining samples containing 'IDH' in the subtype string into an 'LGG_IDH' group.
xcell_score$molecular_subtype_summarized[grep("IDH", xcell_score$molecular_subtype)] <- "LGG_IDH"
# Group all remaining samples containing 'RTK' or 'FGFR' in the subtype string into an 'LGG_RTK' subgroup
xcell_score$molecular_subtype_summarized[grep("RTK|FGFR", xcell_score$molecular_subtype)] <- "LGG_RTK"
# Modify MYB other MAPK
xcell_score$molecular_subtype_summarized[grep("LGG, MYB/MYBL1, other MAPK", xcell_score$molecular_subtype)] <- "LGG, MYB/MYBL1"

# replace old molecular_subtype column with new column
xcell_score <- xcell_score %>%
  dplyr::select(-c(molecular_subtype)) %>%
  dplyr::rename("molecular_subtype" = "molecular_subtype_summarized")

# output xcell score
xcell_score %>%
  write_tsv(file.path(xcell_output_dir, "xcell_score_cluster.tsv"))

# apply v.test function per gene
vtest_output_all <- plyr::ddply(.data = xcell_score, 
                                .variables = "cell_type", 
                                .fun = function(x) compute.v.test(x, clustering_col = "cluster_assigned"))
# output vtest score
vtest_output_all <- vtest_output_all %>%
  readr::write_tsv(file.path(vtest_output_dir, "vtest_scores_all.tsv"))
