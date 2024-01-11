# ssGSEA using the GSVA implementation to generate Reactome pathway scores (curated c2 pathways from msigdb) 
# on the rnaseq data derived from LGG participant list 

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(GSVA)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "lgg_risk_analysis")
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# read histology file 
histology_file <- file.path(data_dir, "20230826_release.annotated_histologies.tsv") %>%
  read_tsv() 

# read imaging risk file and pull corresponding bs identifiers
lgg_clusters <- readxl::read_xlsx(file.path(input_dir, "RiskScores_Grouping.xlsx"))
histology_file <- histology_file %>%
  filter(cohort_participant_id %in% lgg_clusters$SubjectID,
         experimental_strategy == "RNA-Seq") 

# read tpm data
tpm_file <- file.path(data_dir, "20230826_release-gene-expression-rsem-tpm.collapsed.rds")
tpm_mat <- readRDS(tpm_file)
tpm_mat <- tpm_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = file.path(data_dir, "gencode.v39.primary_assembly.annotation.gtf.gz"))
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
human_reactome_gs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") 
human_reactome_gs <- human_reactome_gs %>% 
  dplyr::select(gs_name, human_gene_symbol)
human_reactome_gs <- base::split(human_reactome_gs$human_gene_symbol, list(human_reactome_gs$gs_name))

# ssgsea scores
ssgsea_scores <- GSVA::gsva(expr = as.matrix(log2(tpm_mat + 1)), 
                            gset.idx.list = human_reactome_gs,
                            method = "ssgsea",
                            min.sz = 10, 
                            max.sz = 500,
                            mx.diff = TRUE) ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
ssgsea_scores <- ssgsea_scores %>% as.data.frame()
saveRDS(object = ssgsea_scores, file = file.path(output_dir, "ssgsea_matrix.rds"))
