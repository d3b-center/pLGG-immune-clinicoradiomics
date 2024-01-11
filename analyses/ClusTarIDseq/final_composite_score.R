#' Function to compute final composite score on NB.mclust + CCP chosen clusters
#'
#' @author Komal S. Rathi
#'
#' @import dplyr
#'
#' @param expr_mat path to expression matrix with raw or expected counts
#' @param hist_file path to histology file with survival info
#' @param filter_expr filter genes by expession. logical TRUE or FALSE
#' @param dispersion_percentile_val value for DGCA::filterGenes dispersion_percentile parameter
#' @param protein_coding_only filter genes to protein coding only. logical TRUE or FALSE
#' @param gencode_version gencode version to fetch gtf. Default is 27
#' @param feature_selection feature selection strategy. Either variance or dip.test
#' @param var_prop proportion of most variable genes to choose
#' @param transformation_type transformation method
#' @param ccp_output path to top-ranking result from CCP
#' @param nbmclust_output path to NB.MClust chosen clusters
#' @param dbscan_output path to DBSCAN chosen clusters
#' @param intnmf_output path to IntNMF chosen clusters
#' @param output_dir output directory to write out output
#'
#' @return
#' a tsv file with the following columns:
#' clustering_method, 
#' cox_anova_pval, cox_anova_chisq i.e pvalue and chisq for effect of cluster membership in a cox model and compare to a 'base'/no-covariate model,
#' cluster.size, noisen, diameter, average.distance, median.distance,
#' average.between, average.within, within.cluster.ss, clus.avg.silwidths, avg_sil,
#' dunn, dunn2, entropy, wb.ratio,
#' cluster_qual cluster quality from COINr,
#' rank final rank assigned by COINr
#'
#' @export
#'

final_composite_score <- function(expr_mat, hist_file,
                                  filter_expr = TRUE, dispersion_percentile_val = 0.2,
                                  protein_coding_only = TRUE, gencode_version = 27,
                                  feature_selection = c("variance", "dip.test"),
                                  var_prop = NULL, transformation_type = c("none", "tmm", "vst", "uq", "log2"),
                                  ccp_output, nbmclust_output = NULL, dbscan_output, intnmf_output,
                                  output_dir) {
  # create output directory
  dir.create(output_dir, recursive = T, showWarnings = F)

  # read histology file
  if (!is.null(hist_file)) {
    hist_file <- read.delim(hist_file)
  }

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

  # read clustering output files and combine in a list
  ccp_output <- readr::read_tsv(ccp_output)
  dbscan_output <- readr::read_tsv(dbscan_output)
  intnmf_output <- readr::read_tsv(intnmf_output)
  # nbmclust cannot be applied to non-expression datasets*
  if(!is.null(nbmclust_output)){
    nbmclust_output <- readr::read_tsv(nbmclust_output)
    cluster_output <- list(ccp_output, nbmclust_output, dbscan_output, intnmf_output)
    names(cluster_output) <- c(
      unique(ccp_output$method),
      unique(nbmclust_output$method),
      unique(dbscan_output$method),
      unique(intnmf_output$method)
    )
  } else {
    cluster_output <- list(ccp_output, dbscan_output, intnmf_output)
    names(cluster_output) <- c(
      unique(ccp_output$method),
      unique(dbscan_output$method),
      unique(intnmf_output$method)
    )
  }
  
  # do for each clustering output
  final_output <- data.frame()
  for (i in 1:length(cluster_output)) {
    # get cluster membership
    cluster_class <- cluster_output[[i]]
    max_cluster <- max(cluster_class$cluster_assigned)

    # combine cluster membership with survival data to compute additional stats
    if (!is.null(hist_file) & max_cluster > 1) {
      data_for_stats <- cluster_class %>%
        inner_join(hist_file %>%
          dplyr::select(Kids_First_Biospecimen_ID, OS_status, OS_days), by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
        mutate(OS_status = ifelse(OS_status == "LIVING", 0,
          ifelse(OS_status == "DECEASED", 1, NA)
        ), OS_days = as.numeric(OS_days))

      # compute the effect of cluster membership in a cox model and compare to a 'base'/no-covariate model (via stats::anova)
      results_cox_clust <- survival::coxph(survival::Surv(OS_days, OS_status) ~ cluster_assigned, data = data_for_stats)
      results_cox <- survival::coxph(survival::Surv(OS_days, OS_status) ~ 1, data = data_for_stats)
      cox_anova <- anova(results_cox, results_cox_clust)
      cox_anova_pval <- round(as.numeric(broom::tidy(cox_anova)[2, "p.value"]), digits = 2)
      cox_anova_chisq <- round(as.numeric(broom::tidy(cox_anova)[2, "statistic"]), digits = 2)
    } else {
      cox_anova_pval <- NA
      cox_anova_chisq <- NA
    }

    # compute distance on dip test output and evaluate silhouette coefficients.
    if (max_cluster > 1) {
      expr_mat_dist <- factoextra::get_dist(t(expr_mat), method = "pearson")
      clus_stats <- fpc::cluster.stats(d = expr_mat_dist, clustering = cluster_class$cluster_assigned)
      clus_stats$avg_sil <- mean(clus_stats$clus.avg.silwidths) # calculate mean
      clus_stats <- lapply(clus_stats, FUN = function(x) {
        if (length(x) > 1 & !is.null(x)) {
          x <- toString(round(x, 2))
        } else {
          x <- round(as.numeric(x), digits = 2)
        }
      })
      clus_stats <- t(as.data.frame(unlist(clus_stats)))
      clus_stats <- clus_stats %>%
        as.data.frame() %>%
        dplyr::select(cluster.size, noisen, diameter, average.distance, median.distance, average.between, average.within, within.cluster.ss, clus.avg.silwidths, avg_sil, dunn, dunn2, entropy, wb.ratio)
    } else {
      clus_stats <- data.frame(
        "cluster.size" = NA,
        "noisen" = NA,
        "diameter" = NA,
        "average.distance" = NA,
        "median.distance" = NA,
        "average.between" = NA,
        "average.within" = NA,
        "within.cluster.ss" = NA,
        "clus.avg.silwidths" = NA,
        "avg_sil" = NA,
        "dunn" = NA,
        "dunn2" = NA,
        "entropy" = NA,
        "wb.ratio" = NA
      )
    }

    # combine for all clustering approaches
    tmp_df <- data.frame(cox_anova_pval, cox_anova_chisq, clus_stats, row.names = names(cluster_output)[i])
    final_output <- rbind(final_output, tmp_df)
  }

  # compute final rank for all clustering methods
  # subset columns used for scoring methods
  output_df_sub <- final_output %>%
    dplyr::select(
      average.between, average.within,
      within.cluster.ss, dunn, entropy, avg_sil
    ) %>%
    mutate(avg_sil = ifelse(as.numeric(avg_sil) < 0, 0, avg_sil))
  output_df_sub[] <- sapply(output_df_sub, as.numeric)

  # uCode is required field and uName is optional for coin object
  output_data <- output_df_sub %>%
    tibble::rownames_to_column("uName") %>%
    mutate(uCode = uName)
  
  # Replace 0 values with small values for testing scoring methods
  output_data[output_data == 0] <- 0.000001
  
  # Create metadata table for coin object, assigning direction and weights for the 9-component score
  output_meta <- tibble(
    "Level" = c(rep(1, 6), 2),
    "iName" = c(colnames(output_data)[2:7], "cluster_qual"),
    "iCode" = c(colnames(output_data)[2:7], "cluster_qual"),
    "Parent" = c(rep("cluster_qual", 6), rep(NA, 1)), 
    "Direction" = c(1, -1, -1, 1, -1, 1, 1),
    "Weight" = c(1, 1, 2, .1, .01, 1, 1),
    "Type" = c(rep("Indicator", 6), "Aggregate")
  )
  
  # Assemble COIN object with data and metadata, normalise, aggregate, and get results
  ASEM <- COINr::new_coin(iData = output_data, iMeta = output_meta)
  ASEM <- COINr::Normalise(x = ASEM, dset = "Raw")
  ASEM <- COINr::Aggregate(x = ASEM, dset = "Normalised")
  rslts <- COINr::get_results(coin = ASEM, tab_type = "Summ", dset = "Aggregated")
  rslts <- rslts %>%
    dplyr::rename("rank" = "Rank") # just to be consistent with the general naming convention
  
  # merge with original data and order by ranks
  final_output <- final_output %>% 
    rownames_to_column('clustering_method') %>% 
    inner_join(rslts, by = c("clustering_method" = "uCode")) %>%
    dplyr::arrange(rank) %>%
    dplyr::relocate(clustering_method)

  # save output to file
  dir.create(output_dir, showWarnings = F, recursive = T)
  readr::write_tsv(x = final_output, file = file.path(output_dir, "final_clustering_output.tsv"))
}
