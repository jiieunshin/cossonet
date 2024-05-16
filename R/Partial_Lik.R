#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param time A vector of true response variable.
#' @param status A vector of estimated response variable fitted by `cdcosso`.
#' @param RS A family corresponding to the type of response variable.
#' @param K A family corresponding to the type of response variable.
#' @param a A family corresponding to the type of response variable.
#' @param neg A family corresponding to the type of response variable.
#'
#' @return A vector of Kullback-Leibler divergence evaluated by test data.
#' @export
#'
Partial_Lik = function (time, status, RS, K, a, neg = FALSE) {
  pl = rep(NA, ncol(RS))
  eventtime = unique(time[status == 1])
  tie.size = as.numeric(table(time[status == 1]))
  for (k in 1:ncol(RS)) {
    failid = which(time == eventtime[k])
    pl[k] = tie.size[k] * log(sum(exp(K[RS[,  k],] %*% a), na.rm = T) + 1e-10)
  }
  pl = sum(pl) - t(status) %*% K %*% a

  if(neg) pl = -pl
  return(pl)
}
