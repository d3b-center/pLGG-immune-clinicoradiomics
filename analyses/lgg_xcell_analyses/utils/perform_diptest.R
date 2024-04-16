#' Function: Run diptest and filter matrix based on diptest results
#'
#' @author Run Jin
#' This function runs diptest on all features of a count matrix
#'
#' @import DESeq2
#' @import edgeR
#' @import dplyr
#'
#' @param count_matrix either a raw count matrix or tpm matrix
#' @param matrix_type selected from either "count", "normalized_count" or "tpm"
#' @param min_n mininal number of features to keep. Default value is NULL.
#' if set to NULL, only features passing `pval_cutoff` are kept
#' otherwise, at least `min_n` features with lowest `min_n` pvalues were kept
#' @param pval_cutoff pvalue cutoff for filtering features based on diptest
#' @param normalization_method select from "deseq", "tmm", "vst", "log_cpm", "none"
#' if matrix is "tpm", the normalization method should be "none"
#' @param filter_low_expr TRUE or FALSE parameter to indicate whether lowly expressed genes need to be filtered
#' @param dispersion_percentile_val value for DGCA::filterGenes dispersion_percentile parameter
#' 
#' @return normalized count matrix as a dataframe containing either significant or `min_n` features is returned
#'
#'
#' @export
perform_diptest <- function(count_matrix,
                            matrix_type = "count",
                            min_n = NULL,
                            pval_cutoff = 0.05,
                            normalization_method = c("deseq", "tmm", "vst", "log_cpm", "none"),
                            filter_low_expr = FALSE,
                            dispersion_percentile_val = 0.2) {
  #### add statements to allow only particular combinations
  if (!normalization_method %in% c("deseq", "tmm", "vst", "log_cpm", "none")) {
    stop("Error: Please select a valid normalization method.")
  }
  if (!matrix_type %in% c("tpm", "count", "normalized_count")) {
    stop("Error: Input matrix type is not supported by this function.")
  }
  if (matrix_type == "tpm" && !is.null(normalization_method)) {
    stop("Error: TPM matrix input does not need normalization.")
  }
  if (matrix_type == "count" && normalization_method == "none") {
    stop("Error: Please select a normalization method to normalize the count matrix.")
  }

  ######## filter raw counts if needed
  if (filter_low_expr) {
    count_matrix <- filter_low_expr_gene(count_matrix, dispersion_percentile = dispersion_percentile_val)
  }

  ######## normalize the count
  if (normalization_method == "deseq") {
    # second is DESeq2 normalized count
    dds <- DESeqDataSetFromMatrix(
      countData = round(count_matrix),
      colData = as.data.frame(colnames(count_matrix)),
      design = ~1
    )
    count_matrix_normalized <- estimateSizeFactors(dds)
    count_matrix_normalized <- counts(count_matrix_normalized, normalized = T)
  } else if (normalization_method == "tmm") {
    group <- factor(rep(c(1), times = length(colnames(count_matrix)))) # give them the same group for all genes
    count_matrix_dge <- DGEList(counts = count_matrix, group = group)

    # compute the normalization factors (TMM) for scaling by library size
    count_matrix_normalized <- calcNormFactors(count_matrix_dge, method = "TMM")
    count_matrix_normalized <- cpm(count_matrix_normalized)
  } else if (normalization_method == "vst") {
    count_matrix_normalize <- DESeq2::varianceStabilizingTransformation(round(as.matrix(count_matrix)),
      blind = TRUE,
      fitType = "parametric"
    )
  } else if (normalization_method == "log_cpm") {
    group <- factor(rep(c(1), times = length(colnames(count_matrix)))) # give them the same group for all genes
    count_matrix_dge <- DGEList(counts = count_matrix, group = group)

    # compute the normalization factors (TMM) for scaling by library size
    count_matrix_dge <- calcNormFactors(count_matrix_dge, method = "TMM")
    count_matrix_dge_cpm <- cpm(count_matrix_dge, log = T)
  } else if (normalization_method == "none") {
    count_matrix_normalized <- count_matrix
  }
  
  # define a dataframe to store the pval
  count_matrix_pval <- data.frame()
  count_matrix_normalized <- as.data.frame(count_matrix_normalized)
  for (i in 1:nrow(count_matrix_normalized)) {
    # calculate dip.test for each gene
    diptest_each <- count_matrix_normalized[i, ] %>%
      as.numeric() %>%
      diptest::dip.test(simulate.p.value = FALSE, B = 2000)
    # gather the pvalue for each gene
    pval_each <- diptest_each$p.value
    # add another column to store pval
    count_matrix_pval_each <- count_matrix_normalized[i, ] %>%
      as.data.frame() %>%
      dplyr::mutate(pval = pval_each)
    # cocgine each line back to the dataframe
    count_matrix_pval <- bind_rows(count_matrix_pval, count_matrix_pval_each)
  }

  # define rownames
  rownames(count_matrix_pval) <- rownames(count_matrix_normalized)

  # check how many features pass the p-value threshold
  n <- count_matrix_pval %>% 
    dplyr::filter(pval < pval_cutoff) %>%
    nrow()
  print(n)
  
  if(is.null(min_n)){
    # if min_n is null, only keep significant features 
    count_matrix_pval <- count_matrix_pval %>%
      dplyr::filter(pval < pval_cutoff) %>%
      dplyr::select(-pval)
  } else {
    # if min_n is set to a value, keep at least those number of features
    count_matrix_pval <- count_matrix_pval %>%
      dplyr::arrange(pval, descending = FALSE) %>%
      head(min_n) %>%
      dplyr::select(-pval)
  }
  
  return(count_matrix_pval)
}
