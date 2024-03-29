# Author: Komal S. Rathi
# script to run xCell on samples of interest

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(immunedeconv)
})

# parse arguments
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "input matrix for immune profiling (.rds)"),
  make_option(c("--histology"), type = "character",
              help = "histology file for all OpenPedCan samples (.tsv)"),
  make_option(c("--short_hist"), type = "character", default = NULL,
              help = "short histology filter"),
  make_option(c("--broad_hist"), type = "character", default = NULL,
              help = "broad histology filter"),
  make_option(c("--cg"), type = "character", default = NULL,
              help = "cancer group filter"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
mat <- opt$mat
histology <- opt$histology
short_hist <- opt$short_hist
broad_hist <- opt$broad_hist
cg <- opt$cg
output_dir <- opt$output_dir

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "lgg_xcell_analyses")
dir.create(output_dir, recursive = TRUE, showWarnings = F)

# histology file
histology_df <- readr::read_tsv(histology) %>%
  unique()

# subset to histology of interest
histology_df <- histology_df %>%
  filter(sample_type == "Tumor",
         composition != "Derived Cell Line", 
         !cohort %in% c("TCGA", "GTEx"),
         experimental_strategy == "RNA-Seq",
         !is.na(molecular_subtype)) %>%
  dplyr::filter(short_histology %in% short_hist | 
                  broad_histology %in% broad_hist | 
                  cancer_group %in% cg) 

# read input matrix file and filter
mat <- readRDS(mat)
mat <- mat %>%
  dplyr::select(any_of(histology_df$Kids_First_Biospecimen_ID))
histology_df <- histology_df %>%
  filter(Kids_First_Biospecimen_ID %in% colnames(mat))
print(dim(mat))

# compute xcell scores
deconv_output <- deconvolute(gene_expression = as.matrix(mat), 
                             method = "xcell", 
                             arrays = F)

# remove cell types and convert to matrix
deconv_output <- deconv_output %>%
  dplyr::filter(!cell_type %in% c("Smooth muscle", "Osteoblast", "Myocytes", 
                                  "Platelets", "Adipocytes", "MPP", "Hepatocytes", 
                                  "Skeletal muscle", "Keratinocytes", "Sebocytes", 
                                  "Erythrocytes", "Epithelial cells", "Chondrocytes", 
                                  "ly Endothelial cells", "Endothelial cells", "Neurons")) %>%
  column_to_rownames("cell_type") 

# # convert to long format
# deconv_output <- deconv_output %>%
#   as.data.frame() %>%
#   gather(Kids_First_Biospecimen_ID, fraction, -c(cell_type)) %>%
#   as.data.frame()
# 
# # merge output with clinical data
# full_output <- histology_df %>%
#   inner_join(deconv_output, by = "Kids_First_Biospecimen_ID") %>%
#   mutate(method = "xcell")

# save xcell scores
saveRDS(deconv_output, file = file.path(output_dir, "xcell_scores.rds"))
