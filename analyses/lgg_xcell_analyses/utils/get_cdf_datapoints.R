#' Modified function of ConsensusClusterPlus:::CDF to obtain cdf data points and delta area values
#'
#' @author Komal S. Rathi
#'
#' @param ml list of consensus matrix result for each k
#' @param breaks number of data breaks, default is 100
#'
#' @return a list with a dataframe of cdf data points and a vector of delta auc values
#'
#' @export

get_cdf_datapoints <- function(ml, breaks = 100) {
  
  # set seed for reproducibility
  set.seed(100)
  
  k <- length(ml)
  this_colors <- rainbow(k - 1)
  areaK <- c()
  cdf_data <- list()
  for (i in 2:length(ml)) {
    v <- ConsensusClusterPlus:::triangle(ml[[i]], mode = 1)
    h <- hist(v, plot = FALSE, breaks = seq(0, 1, by = 1 / breaks))
    h$counts <- cumsum(h$counts) / sum(h$counts)
    thisArea <- 0
    for (bi in 1:(length(h$breaks) - 1)) {
      thisArea <- thisArea + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi])
      bi <- bi + 1
    }
    areaK <- c(areaK, thisArea)
    # cdf data points
    cdf_data[[i]] <- h
  }

  # delta AUC
  delta_auc <- areaK[1]
  for (i in 2:(length(areaK))) {
    delta_auc <- c(delta_auc, (areaK[i] - areaK[i - 1]) / areaK[i - 1])
  }

  # return list with CDF data points and delta AUC
  return(list(cdf_data = cdf_data, delta_auc = delta_auc))
}
