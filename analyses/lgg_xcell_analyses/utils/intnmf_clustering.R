#' Function to perform Integrative Clustering of Multiple Genomic Dataset
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
#' @param max_k number of k to evaluate.
#' @param output_dir output directory to write out output
#'
#' @return
#' a tsv file with sample and assigned clusters with the most optimal intNMF clustering
#'
#' @export
#'

intnmf_clustering <- function(expr_mat,
                              filter_expr = TRUE, dispersion_percentile_val = 0.2,
                              protein_coding_only = TRUE, gencode_version = 27,
                              feature_selection = c("variance", "dip.test"),
                              var_prop = NULL, min_n = NULL, transformation_type = c("none", "tmm", "vst", "uq", "log2", "rank"),
                              max_k = NULL, output_dir) {
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

  # do some filtering based on expression
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
  dat <- list(t(expr_mat))
  
  # assign weights to each modality
  wt = if(is.list(dat)) rep(1,length(dat)) else 1
  
  ##### outdated ####
  # # find optimum number of clusters for the data 
  # opt.k <- nmf.opt.k(
  #   dat = dat,
  #   n.runs = 15,
  #   maxiter = 100,
  #   k.range = 2:max_k
  # )
  # saveRDS(opt.k, file = file.path(output_dir, "intnmf_optk.rds"))
  # 
  # # run Nonnegative Matrix Factorization of Multiple data using Nonnegative Alternating Least Square
  # opt_k <- as.numeric((which.max(apply(opt.k, MARGIN = 1, mean))) + 1)
  # nmf_output <- IntNMF::nmf.mnnals(
  #   dat = dat,
  #   k = opt_k,
  #   seed = TRUE,
  #   wt = if (is.list(dat)) rep(1, length(dat)) else 1
  # )
  # saveRDS(nmf_output, file = file.path(output_dir, "intnmf_fit.rds"))
  ##### outdated ####
  
  # run run_clusterstats across all k-values
  # get the nmf output corresponding to the most optimal k
  fname <- file.path(output_dir, "intnmf_best_fit.rds")
  if(!file.exists(fname)){
    nmf_output <- run_clusterstats(dat = dat, 
                                   wt = wt, 
                                   output_dir = output_dir, 
                                   k_value = max_k)
  } else {
    nmf_output = readRDS(fname)
  }

  # ConsensusMatPlot
  # Given the integrative NMF fit object, the function creates image plot of the consensus matrix ordered
  # according to clusters groups. Cleaner block structure indicates stronger clusters
  pdf(file = file.path(output_dir, "intnmf_consensus_plot.pdf"), width = 10, height = 10, onefile = F)
  IntNMF::ConsensusMatPlot(fit = nmf_output, rowLab = TRUE, colLab = TRUE)
  dev.off()

  # SilhouettePlot
  # Silhouette width plot is returned together with mean silhouette width for each group, overall silhouette width and summary statistics.
  pdf(file = file.path(output_dir, "intnmf_silhouette_plot.pdf"), width = 10, height = 10, onefile = F)
  IntNMF::SilhouettePlot(fit = nmf_output, cluster.col = NULL)
  dev.off()

  # output
  df <- data.frame(sample = names(nmf_output$clusters), cluster_assigned = nmf_output$clusters, method = "IntNMF")
  write_tsv(df, file = file.path(output_dir, "intnmf_optimal_clusters.tsv"))
}
