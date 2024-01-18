# Generate heatmap with consensus clustering results + harmonized diagnosis
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(pheatmap)
})

# parse arguments
option_list <- list(
  make_option(c("--xcell_file"), type = "character",
              help = "xcell score + cluster file (.tsv)"),
  make_option(c("--plots_dir"), type = "character",
              help = "output directory to save plots")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
xcell_file <- opt$xcell_file
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# read in xcell score file
xcell_score <- readr::read_tsv(xcell_file)

# molecular subtypes > 5 samples 
mol_subtypes_to_use <- xcell_score %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  unique() %>%
  group_by(molecular_subtype) %>%
  summarise(n = n()) %>%
  filter(n >= 5) %>% 
  pull(molecular_subtype)

# filter to selected molecular subtypes and create a wide format
xcell_spread <- xcell_score %>%
  filter(molecular_subtype %in% mol_subtypes_to_use) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cell_type, fraction) %>%
  tidyr::spread(key = cell_type, value = fraction) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID") %>%
  t() %>%
  as.data.frame()

# arrange the annotation file based on the cluster column 
xcell_score$cluster_assigned <- as.factor(xcell_score$cluster_assigned)
anno_file <- xcell_score %>% 
  filter(molecular_subtype %in% mol_subtypes_to_use) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype, cluster_assigned) %>%
  ungroup() %>%
  distinct() %>% 
  arrange(cluster_assigned, molecular_subtype) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

# arrange the score matrix to match the anno file 
xcell_spread <- xcell_spread %>%
  dplyr::select(rownames(anno_file))

# filter low immune scores
xcell_spread <- xcell_spread %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::filter(row_sum > 0.1) %>%
  dplyr::select(-c("row_sum"))

# read annotation colors
anno_colors <- rlist::list.load(file = file.path(data_dir, "colors.yaml"))
anno_colors <- lapply(anno_colors, function(x) unlist(x))
anno_colors$molecular_subtype <- anno_colors$molecular_subtype[names(anno_colors$molecular_subtype) %in% anno_file$molecular_subtype]

# generate heatmap
pheatmap::pheatmap(as.matrix(xcell_spread),
                   cellheight = 12,
                   annotation_col = anno_file,
                   annotation_colors = anno_colors,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   color = colorRampPalette(c("white", "red"))(100),
                   show_colnames = F,
                   show_rownames = T,
                   fontsize = 10,
                   filename = file.path(plots_dir, "xcell_heatmap.pdf"),
                   width = 12,
                   height = 8)
