#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param y A vector of true response variable.
#' @param fhat A vector of estimated response variable fitted by `cdcosso`.
#' @param obj Argument corresponding to the type of response variable.
#'
#' @return A value of Kullback-Leibler divergence evaluated by tes data.
#' @export
#'

KLD = function(y, fhat, obj){
  if(obj$family == "gaussian") B = function(x) x
  if(obj$family == "binomial") B = function(x) log(exp(x) + 1)
  if(obj$family == "poisson") B = function(x) exp(x)

  return(- y * fhat + B(fhat))
}
