#' Function: Run NB.MClust function and return cluster results for samples
#'
#' @author Run Jin
#' This function runs NB.MClust clustering for features of the raw count matrix
#'
#' @param expr_mat path to expression matrix with raw or expected counts
#' @param filter_expr filter genes by expession. logical TRUE or FALSE
#' @param dispersion_percentile_val value for DGCA::filterGenes dispersion_percentile parameter
#' @param protein_coding_only filter genes to protein coding only. logical TRUE or FALSE
#' @param gencode_version gencode version to fetch gtf. Default is 27
#' @param feature_selection feature selection strategy. Either variance or dip.test
#' @param var_prop proportion of most variable genes to choose
#' @param max_k maximum cluster number to search through
#' @param output_dir output directory to write out output
#'
#' @import dplyr
#'
#' @return a dataframe with one column specifying samples clustered and another one the cluster assigned
#'
#' @export
#'

nbmclust_clustering <- function(expr_mat,
                                filter_expr = TRUE, dispersion_percentile_val = 0.2,
                                protein_coding_only = TRUE, gencode_version = 27,
                                feature_selection = c("variance", "dip.test"),
                                var_prop = NULL,
                                max_k = 10,
                                output_dir) {
  # create output directory
  dir.create(output_dir, showWarnings = F, recursive = T)

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

  # calculate features to use for clustering based on input
  if (feature_selection == "variance") {
    # reduce to the most variable features, measured by median absolute deviation
    var_prop <- as.numeric(var_prop)
    stopifnot(var_prop <= 100)
    print(paste("subset to", var_prop, "%", "variable features"))
    n <- round(var_prop / 100 * nrow(expr_mat))
    mads <- apply(expr_mat, 1, mad)
    expr_mat <- expr_mat[rev(order(mads))[1:n], ]
  } else if (feature_selection == "dip.test") {
    # dip.test to filter features
    expr_mat <- perform_diptest(count_matrix = expr_mat, matrix_type = "normalized_count", normalization_method = "none", filter_low_expr = FALSE)
  }
  print(paste0("Number of features:", nrow(expr_mat)))

  # transpose the matrix for fitting the data
  expr_mat_transformed <- t(expr_mat)

  # round the matrix
  expr_mat <- round(expr_mat)

  # run NB.MClust
  nb_mclust_res <- NB.MClust::NB.MClust(
    Count = expr_mat_transformed,
    K = 2:max_k,
    tau0 = 10,
    rate = 0.9,
    iteration = 50
  )

  # estimated cluster assignment
  cluster_nb <- nb_mclust_res$cluster %>%
    as.data.frame()
  colnames(cluster_nb) <- "cluster_assigned"
  cluster_nb <- cluster_nb %>%
    tibble::rownames_to_column("sample") %>%
    mutate(method = "NB.MClust")

  # save to file
  readr::write_tsv(x = cluster_nb, file = file.path(output_dir, "nbmclust_optimal_clusters.tsv"))
}
