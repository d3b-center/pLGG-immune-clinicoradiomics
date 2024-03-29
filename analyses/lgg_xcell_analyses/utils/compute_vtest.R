#' function to compute v-test for each feature across clusters
#'
#' @author Komal S. Rathi
#'
#' @param x a dataframe containing normalized count of a particular gene (feature) for each cluster
#' @param clustering_col the column name of the column containing the cluster assignment
#'
#' @return a dataframe containing clusters, v-score for the cluster for the particular gene/feature
#' column name for vscore is `v_score` and column name for cluster assignment is `cluster`
#'
#' @export

compute_vtest <- function(x, clustering_col) {
  # get unique clusters i.e. 1-n clusters
  clusters <- unique(x[[clustering_col]])

  out <- data.frame()
  # iterate over each cluster
  for (i in 1:length(clusters)) {
    # subset x to cluster i
    y <- x %>%
      filter(get(clustering_col) == clusters[i])

    # mean feature expression per cluster
    mean.feature.cluster <- unique(y$cluster_feature_mean_score)

    # mean feature expression (i.e. global mean)
    mean.feature <- unique(y$feature_mean_score)

    # calculate numerator
    num <- mean.feature.cluster - mean.feature

    # total sample size
    n <- nrow(x)

    # cluster sample size
    ng <- nrow(y)

    # variance of feature values (i.e. global variance)
    var.feature <- unique(x$feature_variance)

    # calculate denominator
    denom <- (n - ng / n - 1) * (var.feature / ng)
    denom <- sqrt(denom)

    # calculate vscore
    v <- num / denom
    out[i, "cluster"] <- clusters[i]
    out[i, "v_score"] <- v
  }
  return(out)
}
