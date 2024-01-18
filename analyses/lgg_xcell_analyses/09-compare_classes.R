# Function to compare imaging clusters and xCell-derived clusters
suppressPackageStartupMessages({
  library(tidyverse)
  library(mclust)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "lgg_xcell_analyses")
results_dir <- file.path(analysis_dir, "results")

# imaging clusters (imaging data has cohort participant id instead of biospecimen id)
imaging_clusters <- file.path(data_dir, "ImagingClusterAssignment_Aug2023.xlsx") %>%
  readxl::read_xlsx()

# map cohort participant id to biospecimen from histology
histology_df <- file.path(data_dir, "20230826_release.annotated_histologies_subset.tsv") %>% read_tsv()
histology_df <- histology_df %>%
  dplyr::select(Kids_First_Biospecimen_ID, cohort_participant_id) 
imaging_clusters <- imaging_clusters %>%
  dplyr::rename("cohort_participant_id" = "SubjectID",
                "imaging_cluster_assigned" = "Cluster Assignment") %>%
  inner_join(histology_df, by = "cohort_participant_id")

# read optimal output from xCell-derived clustering
optimal_clusters <- readr::read_tsv(file.path(analysis_dir, "results", "xcell_output", "xcell_score_cluster.tsv"))
optimal_clusters <- optimal_clusters %>%
  dplyr::select(Kids_First_Biospecimen_ID, cluster_assigned) %>%
  dplyr::rename("immune_cluster_assigned" = "cluster_assigned") %>%
  unique()

# both classes should have the same biospecimens
df <- imaging_clusters %>%
  inner_join(optimal_clusters)

# rownames are identical
ar_index <- adjustedRandIndex(x = df$imaging_cluster_assigned, y = df$immune_cluster_assigned)
write.table(ar_index, file = file.path(results_dir, "adjusted_rand_index.tsv"),  col.names = F, quote = F, row.names = F)
