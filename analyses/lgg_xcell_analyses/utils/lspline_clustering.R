#' Function to run lspline clustering
#'
#' @author Komal S. Rathi
#'
#' @import dplyr
#'
#' @param expr_mat path to expression matrix with raw or expected counts
#' @param hist_file path to histology file
#' @param filter_expr filter genes by expession. logical TRUE or FALSE
#' @param dispersion_percentile_val value for DGCA::filterGenes dispersion_percentile parameter
#' @param protein_coding_only filter genes to protein coding only. logical TRUE or FALSE
#' @param gencode_version gencode version to fetch gtf. Default is 27
#' @param feature_selection feature selection strategy. Either variance or dip.test
#' @param var_prop proportion of most variable genes to choose
#' @param transformation_type expression matrix transformation strategies.
#' Can be either none for no transformation,
#' vst for variance stabilizing transformation,
#' uq for upper quantile normalization or
#' log2 for log-transformation
#' @param max_k maximum k to evaluate
#' @param coef_cutoff cutoff of data-point threshold to determine longest stretch closest to zero
#' @param min_cluster_size_prop min proportion of samples required in a cluster
#' @param max_cluster_size_prop max proportion of samples allowed in a cluster
#' @param compute_all_equal identify if clusters are equally distributed. logical TRUE or FALSE
#' @param output_dir output directory to write out output for all combinations of inputs
#'
#' @return
#' a tsv file with the following columns:
#' method, algorithm,	distance, k,
#' stretch i.e. length of datapoints closest to 0,
#' delta_auc,	p_val i.e. p-value of ks-test between k and k-1,
#' cox_anova_pval, cox_anova_chisq i.e pvalue and chisq for effect of cluster membership in a cox model and compare to a 'base'/no-covariate model,
#' cluster.size, noisen, diameter, average.distance, median.distance,
#' average.between, average.within, within.cluster.ss, clus.avg.silwidths, avg_sil,
#' dunn, dunn2, entropy, wb.ratio,
#' cluster_size_min TRUE if all clusters have at least the specified prop of samples, FALSE if not or NA when no limit specified,
#' cluster_size_max TRUE if all clusters have at most the specified prop of samples, FALSE if not or NA when no limit specified,
#' all_equal p-value indicating if all clusters are equally distributed,
#' cluster_qual cluster quality from COINr,
#' rank final rank assigned by COINr
#'
#' a tsv file with sample and assigned clusters with the most optimal clustering corresponding to the topmost rank assigned by COINr
#'
#' @export
#'

lspline_clustering <- function(expr_mat, hist_file,
                               filter_expr = TRUE, dispersion_percentile_val = 0.2,
                               protein_coding_only = TRUE, gencode_version = 27,
                               feature_selection = c("variance", "dip.test"),
                               var_prop = NULL, transformation_type = c("none", "tmm", "vst", "uq", "log2"),
                               max_k, coef_cutoff = 0.5,
                               min_cluster_size_prop = NULL, max_cluster_size_prop = NULL,
                               compute_all_equal = TRUE,
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
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(gene %in% gencode_gtf$gene_name) %>%
      tibble::column_to_rownames("gene")
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

  # compute all combinations of input parameters
  algorithms <- c("hc", "pam", "km")
  distances <- c("pearson", "spearman", "euclidean", "manhattan", "binary", "maximum", "canberra", "minkowski")
  compute_all <- expand.grid(algorithms, distances)
  colnames(compute_all) <- c("algorithms", "distances")

  # km can only be used with euclidean
  compute_all <- compute_all %>%
    dplyr::filter(!(!distances == "euclidean" & algorithms == "km"))
  
  # hc and pam are having issues with manhattan and minkowski distances*
  compute_all <- compute_all %>%
    dplyr::filter(!(distances %in% c("minkowski", "manhattan") & algorithms %in% c("pam", "hc")))

  # run ConsensusClusterPlus for all combinations
  set.seed(100)

  # output dataframe
  output_df <- data.frame()
  for (i in 1:nrow(compute_all)) {
    ccp_algorithm <- as.character(compute_all[i, "algorithms"])
    ccp_distance <- as.character(compute_all[i, "distances"])
    print(ccp_algorithm)
    print(ccp_distance)

    # output filename
    output_file <- paste(ccp_algorithm, ccp_distance, var_prop, sep = "_")
    output_file <- paste0(output_file, ".rds")
    output_file <- file.path(output_dir, output_file)

    # call function to cluster via ConsensusClusterPlus
    ccp_output <- run_ccp(
      expr_mat = expr_mat,
      algorithm = ccp_algorithm,
      distance = ccp_distance,
      max_k = max_k,
      output_file = output_file,
      ccp_dir = output_dir
    )

    # extract and create a consensus matrix list for each k
    ml <- list()
    for (j in 2:length(ccp_output)) {
      ml[[j]] <- ccp_output[[j]]$ml
    }

    # get_cdf_datapoints: modified function of ConsensusClusterPlus:::CDF
    # it returns a list of cdf data points and delta area values
    out <- get_cdf_datapoints(ml = ml, breaks = 100)
    cdf_data <- out$cdf_data # cdf data points
    delta_auc <- c(0, out$delta_auc) # delta auc values
    names(delta_auc) <- paste0("k", 1:length(delta_auc))

    # function to compute longest stretch closest to 0
    longestNAstrech <- function(x) {
      with(rle(is.na(x)), max(lengths[values]))
    }

    # compute the longest stretch for each k for the selected distance + algorithm combination
    n <- length(delta_auc)
    coefs <- list()
    k_flat <- data.frame() # output dataframe
    for (item in 2:n) {
      k <- as.numeric(gsub("k", "", names(delta_auc)[item]))
      # taken from ConsensusClusterPlus:::CDF, plot original consensus CDF plot
      x <- cdf_data[[k]]$mids
      y <- cdf_data[[k]]$counts
      m1 <- lm(y ~ lspline::lspline(x, seq(0.005, 0.995, 0.02)))

      # get coefficients that are below coef_cutoff and convert to NA to calculate longest stretch
      coef <- m1$coefficients
      coef[coef <= coef_cutoff] <- NA

      # get cluster membership for a given k
      cluster_class <- ccp_output[[item]]$consensusClass
      cluster_class <- data.frame(cluster_class)

      # combine cluster membership with survival data to compute additional stats
      if (!is.null(hist_file)) {
        data_for_stats <- cluster_class %>%
          tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
          inner_join(hist_file %>%
            dplyr::select(Kids_First_Biospecimen_ID, OS_status, OS_days), by = "Kids_First_Biospecimen_ID") %>%
          mutate(OS_status = ifelse(OS_status == "LIVING", 0,
            ifelse(OS_status == "DECEASED", 1, NA)),
            OS_days = as.numeric(OS_days))

        # compute the effect of cluster membership in a cox model and compare to a 'base'/no-covariate model (via stats::anova)
        results_cox_clust <- survival::coxph(survival::Surv(OS_days, OS_status) ~ cluster_class, data = data_for_stats)
        results_cox <- survival::coxph(survival::Surv(OS_days, OS_status) ~ 1, data = data_for_stats)
        cox_anova <- anova(results_cox, results_cox_clust)
        cox_anova_pval <- round(as.numeric(broom::tidy(cox_anova)[2, "p.value"]), digits = 2)
        cox_anova_chisq <- round(as.numeric(broom::tidy(cox_anova)[2, "statistic"]), digits = 2)
      } else {
        cox_anova_pval <- NA
        cox_anova_chisq <- NA
      }

      # compute distance on dip test output and evaluate silhouette coefficients.
      expr_mat_dist <- factoextra::get_dist(t(expr_mat), method = "pearson")
      clus_stats <- fpc::cluster.stats(d = expr_mat_dist, cluster_class$cluster_class)
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

      # get cluster membership count for a given k
      cluster_class <- cluster_class %>%
        group_by(cluster_class) %>%
        summarise(cluster_size = n())

      # limits on cluster size
      if (!is.null(min_cluster_size_prop) & !is.null(max_cluster_size_prop)) {
        # set limits on cluster size
        # max number of samples allowed in each cluster
        max_cluster_size <- max_cluster_size_prop * ncol(expr_mat)
        # min number of samples required in each cluster
        min_cluster_size <- min_cluster_size_prop * ncol(expr_mat)
        # return TRUE if ALL clusters in the selected k follow the above rule, else return FALSE
        cluster_size <- cluster_class %>%
          mutate(
            cluster_size_min = ifelse(cluster_size >= min_cluster_size, TRUE, FALSE),
            cluster_size_max = ifelse(cluster_size <= max_cluster_size, TRUE, FALSE)
          ) %>%
          summarise(
            cluster_size_min = all(cluster_size_min),
            cluster_size_max = all(cluster_size_max)
          ) %>%
          as.data.frame()
      } else {
        # no limits on cluster size
        cluster_size <- data.frame(cluster_size_min = NA, cluster_size_max = NA)
      }

      # chisq test on sample count per cluster within a k to determine equality in cluster membership
      if (compute_all_equal) {
        all_equal <- chisq.test(x = cluster_class$cluster_size)
        all_equal <- format(all_equal$p.value, digits = 3, scientific = TRUE)
      } else {
        all_equal <- NA
      }

      # calculate ks-test p-value between adjacent k and k-1
      # should conflict arise, this will tell if there is a significant difference between 2 adjacent k
      if (k == 2) {
        # there is no k = 1 so use 0 as baseline
        pval <- ks.test(x = 0, y = cdf_data[[k]]$counts, alternative = "two.sided")
      } else {
        pval <- ks.test(
          x = cdf_data[k - 1][[1]]$counts,
          y = cdf_data[[k]]$counts,
          alternative = "two.sided"
        )
      }

      # create a dataframe of all k, longest stretch of data points close to 0, delta AUC and p-value between k and k-1
      # cluster_size and all_equal filters
      k_df <- data.frame(
        k = k,
        stretch = longestNAstrech(coef),
        delta_auc = round(delta_auc[[item]], digits = 2),
        p_val = format(pval$p.value, digits = 3, scientific = TRUE),
        cox_anova_pval, cox_anova_chisq,
        clus_stats,
        cluster_size,
        all_equal
      )

      # bind for all k
      k_flat <- rbind(k_flat, k_df)
    }

    # assign the output to corresponding algorithm and distance
    tmp_df <- data.frame(algorithm = ccp_algorithm, distance = ccp_distance, k_flat)

    # this dataframe with contain output for all combinations of distance + algorithm + k
    output_df <- rbind(output_df, tmp_df)
  }

  # assign rownames
  rownames(output_df) <- paste(output_df$algorithm, output_df$distance, output_df$k, sep = "_")

  # subset columns used for scoring methods
  output_df_sub <- output_df %>%
    dplyr::select(
      stretch, delta_auc, p_val, average.between, average.within,
      within.cluster.ss, dunn, entropy, avg_sil
    ) %>%
    mutate(avg_sil = ifelse(as.numeric(avg_sil) < 0, 0, avg_sil))
  output_df_sub[] <- sapply(output_df_sub, as.numeric)

  # COINr results
  
  # uCode is required field and uName is optional for coin object
  output_data <- output_df_sub %>%
    tibble::rownames_to_column("uName") %>%
    mutate(uCode = uName)

  # Replace 0 values with small values for testing scoring methods
  output_data[output_data == 0] <- 0.000001

  # Create metadata table for coin object, assigning direction and weights for the 9-component score
  output_meta <- tibble(
  "Level" = c(rep(1, 9), 2),
  "iName" = c(colnames(output_data)[2:10], "cluster_qual"),
  "iCode" = c(colnames(output_data)[2:10], "cluster_qual"),
  "Parent" = c(rep("cluster_qual", 9), rep(NA, 1)), 
  "Direction" = c(1, 1, -1, 1, -1, -1, 1, -1, 1, 1),
  "Weight" = c(3, 1, 1, 1, 1, 2, .1, .01, 1, 1),
  "Type" = c(rep("Indicator", 9), "Aggregate")
  )

  # Assemble COIN object with data and metadata, normalise, aggregate, and get results
  ASEM <- COINr::new_coin(iData = output_data, iMeta = output_meta)
  ASEM <- COINr::Normalise(x = ASEM, dset = "Raw")
  ASEM <- COINr::Aggregate(x = ASEM, dset = "Normalised")
  rslts <- COINr::get_results(coin = ASEM, tab_type = "Summ", dset = "Aggregated")
  rslts <- rslts %>%
    dplyr::rename("rank" = "Rank") # just to be consistent with the general naming convention

  # add to original dataframe and arrange by coinr_rank
  output_df <- output_df %>%
    tibble::rownames_to_column("method") %>%
    inner_join(rslts, by = c("method" = "uCode")) %>%
    arrange(rank)

  # write filtered output file to be evaluated
  readr::write_tsv(x = output_df, file = file.path(output_dir, "lspline_output.tsv"))

  # pick the top ranking method and write out clustering info
  output_file <- paste(output_df[1, "algorithm"], output_df[1, "distance"], var_prop, sep = "_")
  output_file <- paste0(output_file, ".rds")
  output_file <- file.path(output_dir, output_file)
  optimal_clusters <- readRDS(output_file)
  optimal_clusters <- optimal_clusters[[output_df[1, "k"]]]$consensusClass
  optimal_clusters <- data.frame(cluster_assigned = optimal_clusters)
  optimal_clusters <- optimal_clusters %>%
    tibble::rownames_to_column("sample") %>%
    mutate(method = output_df[1, "method"])
  readr::write_tsv(x = optimal_clusters, file = file.path(output_dir, "ccp_optimal_clusters.tsv"))
}
