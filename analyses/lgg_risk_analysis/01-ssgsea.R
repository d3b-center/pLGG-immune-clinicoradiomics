# ssGSEA using the GSVA implementation to generate Reactome pathway scores (curated c2 pathways from msigdb)
# on the rnaseq data derived from LGG participant list

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(msigdbr)
  library(GSVA)
})

# parse arguments
option_list <- list(
  make_option(c("--mat"), type = "character", help = "input RNA-seq gene expression matrix (preferably TPM)"),
  make_option(c("--histology_file"), type = "character", help = "histology file for all OpenPedCan samples (.tsv)"),
  make_option(
    c("--risk_file"),
    type = "character",
    default = NULL,
    help = "file with risk scores"
  ),
  make_option(
    c("--gtf_file"),
    type = "character",
    default = NULL,
    help = "gencode GTF file"
  ),
  make_option(c("--output_dir"), type = "character", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
mat <- opt$mat
histology_file <- opt$histology_file
risk_file <- opt$risk_file
gtf_file <- opt$gtf_file
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# read histology file
histology_file <- read_tsv(histology_file)

# read imaging risk file and pull corresponding bs identifiers
lgg_clusters <- readxl::read_xlsx(risk_file)
histology_file <- histology_file %>%
  filter(
    cohort_participant_id %in% lgg_clusters$SubjectID,
    experimental_strategy == "RNA-Seq"
  )

# read tpm data
tpm_mat <- readRDS(mat)
tpm_mat <- tpm_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))

# read gtf and filter to protein coding
gencode_gtf <- rtracklayer::import(con = gtf_file)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()
tpm_mat <- tpm_mat %>%
  rownames_to_column('gene') %>%
  filter(gene %in% gencode_gtf$gene_name) %>%
  column_to_rownames('gene')

# reactome gene set
human_reactome_gs <- msigdbr::msigdbr(species = "Homo sapiens",
                                      category = "C2",
                                      subcategory = "CP:REACTOME")
human_reactome_gs <- human_reactome_gs %>%
  dplyr::select(gs_name, human_gene_symbol)
human_reactome_gs <- base::split(human_reactome_gs$human_gene_symbol,
                                 list(human_reactome_gs$gs_name))

# ssgsea scores
ssgsea_scores <- GSVA::gsva(
  expr = as.matrix(log2(tpm_mat + 1)),
  gset.idx.list = human_reactome_gs,
  method = "ssgsea",
  min.sz = 10,
  max.sz = 500,
  mx.diff = TRUE
) ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
ssgsea_scores <- ssgsea_scores %>% as.data.frame()
saveRDS(object = ssgsea_scores,
        file = file.path(output_dir, "ssgsea_matrix.rds"))
