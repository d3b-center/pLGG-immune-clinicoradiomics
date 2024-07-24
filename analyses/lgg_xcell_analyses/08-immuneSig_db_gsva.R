# GSVA Analysis for lgg cluster data

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(msigdbr)
  library(ggplot2)
  library(GSEABase)
  library(GSVA)
  library(limma)
  library(gplots)
})

# parse arguments
option_list <- list(
  make_option(c("--mat"), type = "character", help = "gene expression matrix preferably TPM (.rds)"),
  make_option(c("--xcell_file"), type = "character", help = "xcell score + cluster file (.tsv)"),
  make_option(c("--output_dir"), type = "character", help = "output directory for files"),
  make_option(c("--plots_dir"), type = "character", help = "output directory for plots")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
mat <- opt$mat
xcell_file <- opt$xcell_file
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# Imaging clusters (imaging data has cohort participant id instead of biospecimen id)
lgg_clusters <-
  file.path(xcell_file) %>%
  readr::read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, cluster_assigned)

lgg_clusters <- lgg_clusters %>%
  distinct(Kids_First_Biospecimen_ID, cluster_assigned, .keep_all = TRUE) %>%
  column_to_rownames(var = "Kids_First_Biospecimen_ID")

# Expression data use @KR script to fix the gene symbols and duplicates
rna_seq_tpm <-
  file.path(mat) %>%
  readRDS()

# subset samples present in cluster data
rna_seq_tpm <- rna_seq_tpm %>%
  dplyr::select(rownames(lgg_clusters))

# Create an expression set
phenoData <- new("AnnotatedDataFrame", data = lgg_clusters)
eset <- Biobase::ExpressionSet(assayData = as.matrix(rna_seq_tpm), phenoData = phenoData)

# Use getGmt() function to read the file
immunesigdb <- msigdbr(category = "C7", subcategory = "IMMUNESIGDB")
immunesigdb <- immunesigdb %>% dplyr::select(gs_name, human_gene_symbol)
immunesigdb <- base::split(immunesigdb$human_gene_symbol, list(immunesigdb$gs_name))

# GSVA analysis (took ~5hrs to calculate absolute values from ranks)
gsva_result <-
  GSVA::gsva(
    eset,
    immunesigdb,
    method = "ssgsea",
    min.sz = 20,
    max.sz = 500
  )
save(gsva_result,
     file = file.path(output_dir, "gsva_igg_cluster_result.RData"))

# Differential expression analysis
design <-
  model.matrix(~ factor(pData(gsva_result)$cluster_assigned))
colnames(design) <- c("cluster1", "cluster2", "cluster3")
fit <- limma::lmFit(exprs(gsva_result), design)
fit <- limma::eBayes(fit)
res <- limma::decideTests(fit, p.value = 0.01)
write.table(res,
            file = file.path(output_dir, "limma_test_results.txt"),
            sep = "\t")

# Create a heatmap for the top 50 differentially expressed pathways
tt <- topTable(fit,
               number = 50,
               adjust.method = "BH",
               sort.by = "B")
DEpwys <- rownames(tt)
DEpwys_es <- exprs(gsva_result[DEpwys, ])

lgg_clusters <-
  lgg_clusters %>% arrange(as.factor(cluster_assigned))

lgg_clusters <- lgg_clusters %>%
  mutate(
    cluster_assigned = case_when(
      cluster_assigned == 1 ~ "cluster1",
      cluster_assigned == 2 ~ "cluster2",
      cluster_assigned == 3 ~ "cluster3",
      TRUE ~ as.character(cluster_assigned)
    )
  )


DEpwys_es <-
  DEpwys_es[, match(rownames(lgg_clusters), colnames(DEpwys_es))]

# Output heatmap
pheatmap::pheatmap(
  DEpwys_es,
  annotation_col = lgg_clusters,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = bluered(256),
  width = 15,
  height = 8,
  show_colnames = FALSE,
  show_rownames = TRUE,
  legend = FALSE,
  scale = 'row',
  main = "Top 50 differentially expressed pathways",
  filename = file.path(plots_dir, "gsva_heatmap.pdf")
)
