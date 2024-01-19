#' Function to perform dbSCAN clustering
#'
#' @author Komal S. Rathi
#'
#' @import dplyr
#' @import kneedle
#'
#' @param expr_mat path to expression matrix with raw or expected counts
#' @param filter_expr filter genes by expession. logical TRUE or FALSE
#' @param dispersion_percentile_val value for DGCA::filterGenes dispersion_percentile parameter
#' @param protein_coding_only filter genes to protein coding only. logical TRUE or FALSE
#' @param gencode_version gencode version to fetch gtf. Default is 27
#' @param feature_selection feature selection strategy. Either variance or dip.test
#' @param var_prop proportion of most variable genes to choose
#' @param transformation_type transformation method
#' @param minpts_val number of minimum points required in the eps neighborhood.
#' @param output_dir output directory to write out output
#'
#' @return
#' a tsv file with sample and assigned clusters with the most optimal dbscan clustering
#'
#' @export
#'

dbscan_clustering <- function(expr_mat,
                              filter_expr = TRUE, dispersion_percentile_val = 0.2,
                              protein_coding_only = TRUE, gencode_version = 27,
                              feature_selection = c("variance", "dip.test"),
                              var_prop = NULL, transformation_type = c("none", "tmm", "vst", "uq", "log2"),
                              minpts_val = NULL, output_dir) {
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
    expr_mat <- DESeq2::varianceStabilizingTransformation(round(as.matrix(expr_mat + 1)), blind = TRUE, fitType = "parametric")
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
    mads <- apply(expr_mat, 1, mad)
    expr_mat <- expr_mat[rev(order(mads))[1:n], ]
  } else if (feature_selection == "dip.test") {
    # dip.test to filter features
    expr_mat <- perform_diptest(count_matrix = expr_mat, matrix_type = "normalized_count", normalization_method = "none", filter_low_expr = FALSE)
  }
  print(paste0("Number of features:", nrow(expr_mat)))

  # get eps value
  x <- 1:ncol(expr_mat)
  y <- sort(dbscan::kNNdist(x = t(expr_mat), k = round(ncol(expr_mat) / 10)))
  knee <- kneedle::kneedle(x, y, sensitivity = 10) # use kneedle algorithm to identify knee/elbow*
  eps <- knee[2] # pull y value for epsilon parameter

  # get minPts value
  if (is.null(minpts_val)) {
    minpts_val <- round(ncol(expr_mat) / 10)
  }

  # run dbscan
  dbscan_out <- dbscan::dbscan(x = t(expr_mat), eps = eps, minPts = minpts_val)
  print(dbscan_out)

  # epsilon plot
  p <- factoextra::fviz_cluster(dbscan_out,
    data = t(expr_mat), stand = FALSE,
    ellipse = TRUE, show.clust.cent = FALSE,
    geom = "point", palette = "jco", ggtheme = ggplot2::theme_classic()
  )
  ggplot2::ggsave(filename = file.path(output_dir, "cluster_plot.pdf"), plot = p)


  # elbow plot
  pdf(file = file.path(output_dir, "KNNdistplot_output.pdf"))
  dbscan::kNNdistplot(t(expr_mat), k = round(ncol(expr_mat) / 10))
  dev.off()

  # output sample to cluster mapping
  optimal_clusters <- data.frame(
    sample = colnames(expr_mat),
    cluster_assigned = dbscan_out$cluster,
    method = "dbscan"
  )
  readr::write_tsv(x = optimal_clusters, file = file.path(output_dir, "dbscan_optimal_clusters.tsv"))
}
