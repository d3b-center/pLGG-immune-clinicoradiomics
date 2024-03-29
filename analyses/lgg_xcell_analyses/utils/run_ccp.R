#' Function to run ConsensusClusterPlus
#'
#' @author Komal S. Rathi
#'
#' @import ConsensusClusterPlus
#'
#' @param expr_mat expression matrix with transformed counts
#' @param algorithm any of the algorithms accepted by ConsensusClusterPlus function for e.g.
#' hc for heirarchical,
#' pam for partition across mediods or
#' km for k-means algorithms
#' @param distance any of the distances accepted by ConsensusClusterPlus function for e.g.
#' pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski or
#' a custom distance function
#' @param max_k maximum k to evaluate
#' @param output_file output filename for the full output of ConsensusClusterPlus to be written out
#' @param ccp_dir output directory of the consensus.pdf output file
#'
#' @return the consensus pdf file and an rds object with full clustering information
#'
#' @export

run_ccp <- function(expr_mat, algorithm = c("hc", "pam", "km"), distance = c("pearson", "spearman", "euclidean", "manhattan", "binary", "maximum", "canberra", "minkowski"), max_k, output_file, ccp_dir) {
  # create directory
  dir.create(ccp_dir, showWarnings = F, recursive = T)

  # set seed for reproducibility
  set.seed(100)
  
  # run ccp
  if(!file.exists(output_file)){
    res_ccp <- invisible(ConsensusClusterPlus::ConsensusClusterPlus(
      d = as.matrix(expr_mat),
      clusterAlg = algorithm,
      finalLinkage = "average",
      distance = distance,
      plot = "pdf",
      writeTable = FALSE,
      reps = 50,
      maxK = max_k,
      pItem = 0.8,
      title = file.path(ccp_dir),
      seed = 100
    ))
    
    # rename consensus output
    new_name <- gsub(".rds", ".pdf", output_file)
    old_name <- file.path(ccp_dir, "consensus.pdf")
    file.rename(from = old_name, to = new_name)
    
    # save ccp output
    saveRDS(res_ccp, file = output_file)
  } else {
    print("File exists")
    res_ccp <- readRDS(file = output_file)
  }
  
  return(res_ccp)
}
