suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(dplyr)
  library(data.table)
  library(MASS)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "lgg_risk_analysis")
input_dir <- file.path(analysis_dir, "input")

# parse arguments
option_list <- list(
  make_option(c("--histology_file"), type = "character", help = "histology file for all OpenPedCan samples (.tsv)"),
  make_option(
    c("--risk_file"),
    type = "character",
    default = NULL,
    help = "file with risk scores"
  ),
  make_option(
    c("--germline_file"),
    type = "character",
    default = NULL,
    help = "file with germline data"
  ),
  make_option(c("--output_dir"), type = "character", help = "output directory"),
  make_option(c("--plots_dir"), type = "character", help = "plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
output_dir <- opt$output_dir
plots_dir <- opt$plots_dir

# output directory
dir.create(output_dir, showWarnings = F, recursive = T)
dir.create(plots_dir, showWarnings = F, recursive = T)

# read histology file
histology_file <- file.path(opt$histology_file) %>%
  fread()

# read risk scores
lgg_risk_groups <- readxl::read_xlsx(opt$risk_file)

# read germline data
lgg_germline <- file.path(opt$germline_file) %>%
  fread() %>%
  dplyr::select(!entrez)

# remove genes with no germline variants
lgg_germline <- lgg_germline %>%
  tibble::column_to_rownames('gene_symbol')

lgg_germline[lgg_germline > 1] <- 1
lgg_germline$total <- rowSums(lgg_germline)

lgg_germline <- lgg_germline %>%
  dplyr::filter(total >= 10) %>%
  dplyr::select(!total) %>%
  t() %>%
  as.data.table(keep.rownames = TRUE) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = rn)

# pull bs_ids for lgg_risk_groups
histology_file <- histology_file %>%
  dplyr::filter(cohort_participant_id %in% lgg_risk_groups$SubjectID,
                sample_type == 'Normal')

lgg_risk_groups <- lgg_risk_groups %>%
  dplyr::inner_join(dplyr::select(histology_file, any_of(
    c('cohort_participant_id', 'Kids_First_Biospecimen_ID')
  )),
  by = c("SubjectID" = "cohort_participant_id"))

# merge risk groups with lgg_germline
lgg_germline <- lgg_germline %>%
  dplyr::inner_join(dplyr::select(lgg_risk_groups, any_of(
    c('Risk Group', 'Kids_First_Biospecimen_ID')
  )),
  by = c("Kids_First_Biospecimen_ID"))

lgg_germline <- lgg_germline %>%
  dplyr::rename(risk_group = `Risk Group`) %>%
  dplyr::mutate(risk_group = as.factor(risk_group)) %>%
  dplyr::mutate(risk_group = factor(
    risk_group,
    levels = c('Low', 'Medium', 'High'),
    ordered = FALSE
  )) %>%
  dplyr::mutate(risk_group_two = ifelse(risk_group == 'Low', 'Low', 'Medium-High'))

genes_of_interest = colnames(lgg_germline)[!colnames(lgg_germline) %in% c('Kids_First_Biospecimen_ID', 'risk_group', 'risk_group_two')]
result = list()

for (i in 1:length(genes_of_interest)) {
  lgg_germline_counts <- lgg_germline %>%
    dplyr::group_by(risk_group, get(genes_of_interest[i])) %>%
    dplyr::summarise(y = n())
  
  fish_res <- fisher.test(
    x = xtabs(
      lgg_germline_counts$y ~ lgg_germline_counts$risk_group + lgg_germline_counts$`get(genes_of_interest[i])`
    ),
    hybrid = T
  )
  result[[i]] <- fish_res
  
}
names(result) <- genes_of_interest

saveRDS(result, file.path(output_dir, 'germline_fisher_test.rds'))
