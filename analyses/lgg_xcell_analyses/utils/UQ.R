#' Function for upper quartile (UQ) normalization function
#'
#' @author Komal S. Rathi
#' 
#' @param X expression matrix  
#'
#' @return normalized matrix
#' @export
#'
#' @examples
#' UQ(expr_mat)
UQ <- function(X){
  
  uq <- function(y){
    quantile(y, 0.75)
  }
  X <- X + 0.1
  upperQ <- apply(X, 2, uq)
  f <- upperQ/mean(upperQ) # calculate normalization factor
  res <- scale(X, center = FALSE, scale = f) 
  return(res)
}