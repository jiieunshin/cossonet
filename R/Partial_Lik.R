#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param time A vector of obbseved time.
#' @param status A vector of observed event indicator.
#' @param K Kernel function evaluated in every component.
#' @param a A vector evaluating observation.
#'
#' @return Partial log-likelihood of Cox proportional model.
#' @export
#'
Partial_Lik = function (time, status, K, a) {
  RS = RiskSet(time, status)
  pl = rep(NA, ncol(RS))
  eventtime = unique(time[status == 1])
  tie.size = as.numeric(table(time[status == 1]))
  for (k in 1:ncol(RS)) {
    failid = which(time == eventtime[k])
    pl[k] = tie.size[k] * log(sum(exp(K[RS[,  k],] %*% a), na.rm = T) + 1e-10)
  }
  pl = sum(pl) - t(status) %*% K %*% a

  return(pl)
}

