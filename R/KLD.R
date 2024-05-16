#' Compute Kullback-Leibler Divergence
#'
#' @param y A vector of true response variable.
#' @param fhat A vector of estimated response variable fitted by `cdcosso`.
#'
#' @return A vector of Kullback-Leibler divergence evaluated by test data.
#'
KLD = function(y, fhat){
  return(- y * fhat + exp(fhat))
}
