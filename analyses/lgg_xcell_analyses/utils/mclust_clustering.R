#' Function to perform mclust clustering
#'
#' @author Komal S. Rathi
#'
#' @import dplyr
#'
#' @param expr_mat path to expression matrix with raw or expected counts
#' @param filter_expr filter genes by expession. logical TRUE or FALSE
#' @param dispersion_percentile_val value for DGCA::filterGenes dispersion_percentile parameter
#' @param protein_coding_only filter genes to protein coding only. logical TRUE or FALSE
#' @param gencode_version gencode version to fetch gtf. Default is 27
#' @param feature_selection feature selection strategy. Either variance or dip.test
#' @param var_prop proportion of most variable genes to choose
#' @param min_n minimum number of genes for dip.test
#' @param transformation_type transformation method
#' @param output_dir output directory to write out output
#'
#' @return
#' a tsv file with sample and assigned clusters with the most optimal mclust clustering
#'
#' @export
#'

mclust_clustering <- function(expr_mat,
                              filter_expr = TRUE, dispersion_percentile_val = 0.2,
                              protein_coding_only = TRUE, gencode_version = 27,
                              feature_selection = c("variance", "dip.test"),
                              var_prop = NULL, min_n = NULL, transformation_type = c("none", "tmm", "vst", "uq", "log2", "rank"),
                              output_dir){
  
  # create output directory
  dir.create(output_dir, recursive = T, showWarnings = F)
  
  # expression data (raw/expected counts)
  expr_mat <- readRDS(expr_mat)
  
  # make sure it is a data-frame
  expr_mat <- as.data.frame(expr_mat)
  
  # remove rows where all counts are 0
  expr_mat <- expr_mat[rowSums(expr_mat[]) > 0, ]
  
  # remove columns and rows with zero standard dev (usually this creates an issue for clustering)
  expr_mat <- expr_mat[apply(expr_mat, MARGIN = 1, function(x) {
    sd(x) != 0
  }), ]
  expr_mat <- expr_mat[, apply(expr_mat, MARGIN = 2, function(x) {
    sd(x) != 0
  })]
  
  # filter based on expression
  if (filter_expr) {
    print("filter by expression")
    expr_mat <- DGCA::filterGenes(
      inputMat = expr_mat,
      filterTypes = c("central", "dispersion"),
      filterDispersionType = "cv",
      filterDispersionPercentile = dispersion_percentile_val,
      sequential = TRUE
    )
  }
  
  # filter to protein coding genes only
  if (protein_coding_only) {
    print("filter to protein coding genes")
    fname <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_", gencode_version, "/gencode.v", gencode_version, ".primary_assembly.annotation.gtf.gz")
    gencode_gtf <- rtracklayer::import(con = fname)
    gencode_gtf <- as.data.frame(gencode_gtf)
    gencode_gtf <- gencode_gtf %>%
      dplyr::select(gene_id, gene_name, gene_type) %>%
      dplyr::filter(gene_type == "protein_coding") %>%
      unique()
    expr_mat <- expr_mat %>%
      rownames_to_column("gene") %>%
      dplyr::filter(gene %in% gencode_gtf$gene_name) %>%
      column_to_rownames("gene")
  }
  
  # transform input matrix
  if (transformation_type == "vst") {
    # transform using vst
    expr_mat <- DESeq2::varianceStabilizingTransformation(round(as.matrix(expr_mat)), blind = TRUE, fitType = "parametric")
    expr_mat <- as.data.frame(expr_mat)
  } else if (transformation_type == "log2") {
    # log-transform
    expr_mat <- log2(expr_mat + 1)
  } else if (transformation_type == "uq") {
    # upper quartile (UQ) normalization
    expr_mat <- UQ(X = expr_mat)
    expr_mat <- as.data.frame(expr_mat)
  } else if (transformation_type == "tmm") {
    # build DGEList and TMM normalize
    y <- edgeR::DGEList(counts = expr_mat)
    y <- edgeR::calcNormFactors(object = y, method = "TMM")
    expr_mat <- as.data.frame(edgeR::cpm(y))
  } else if(transformation_type == "rank") {
    expr_mat <- datawizard::ranktransform(expr_mat %>% as.data.frame())
  } else if (transformation_type == "none") {
    # do nothing
  }
  
  # calculate features to use for clustering based on input
  if (feature_selection == "variance") {
    # reduce to the most variable features, measured by median absolute deviation
    var_prop <- as.numeric(var_prop)
    stopifnot(var_prop <= 100)
    print(paste("subset to", var_prop, "%", "variable features"))
    n <- round(var_prop / 100 * nrow(expr_mat))
    vars <- apply(expr_mat, 1, var)
    expr_mat <- expr_mat[rev(order(vars))[1:n], ]
  } else if (feature_selection == "dip.test") {
    # dip.test to filter features
    expr_mat <- perform_diptest(count_matrix = expr_mat, 
                                matrix_type = "normalized_count", 
                                normalization_method = "none", 
                                filter_low_expr = FALSE, 
                                min_n = min_n)
  }
  print(paste0("Number of features:", nrow(expr_mat)))

  # compute BIC
  BIC <- mclust::mclustBIC(data = t(expr_mat))
  pdf(file = file.path(output_dir, "mclust_bic_plot.pdf"))
  mclust::plot.mclustBIC(x = BIC)
  dev.off()
  
  # run Mclust on computed BIC
  mod1 <- mclust::Mclust(data = t(expr_mat), x = BIC)
  
  # generate and save summary
  summ <- mclust::summary.Mclust(object = mod1, parameters = TRUE)
  saveRDS(summ, file = file.path(output_dir, "mclust_summary.rds"))
  
  # # add classification plot - Error in plot.new() : figure margins too large
  # pdf(file = file.path(output_dir, "mclust_classification.pdf"))
  # mclust::plot.Mclust(mod1, what = "classification")
  # dev.off()
  
  # integrated complete likelihood (summary statistics)
  ICL <- mclust::mclustICL(data = t(expr_mat))
  summ_icl <- mclust::summary.mclustICL(object = ICL)
  saveRDS(summ_icl, file = file.path(output_dir, "mclust_icl_summary.rds"))
  
  # integrated complete likelihood (plot)
  pdf(file = file.path(output_dir, "mclust_icl_plot.pdf"))
  mclust::plot.mclustICL(x = ICL)
  dev.off()
  
  # summary of the likelihood ratio test - takes forever to run even with nboot set to 100 instead of 999
 #  LRT <- mclust::mclustBootstrapLRT(data = t(expr_mat), modelName = "VEI", nboot = 100)
  # saveRDS(LRT, file = file.path(output_dir, "mclust_lrt.rds"))
  
  # output clusters 
  cluster_output <- mod1$classification %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    dplyr::rename("cluster_assigned" = ".") %>%
    mutate(method = "mclust")
  readr::write_tsv(cluster_output, file = file.path(output_dir, "mclust_optimal_clusters.tsv"))
}
