#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param true Explanation variable matrix or data frame.
#' @param est Dependent variable vector or matrix or data frame containing time and status columns (for Cox model).
#'
#' @return A list containing information about the fitted model. Depending on the type of dependent variable, various information may be returned.
#' @export
#'
metric = function(true, est){
  result_tab = table(true, est)
  if(dim(result_tab)[2] == 1){
    result_tab = cbind(0, result_tab)}

  # result.tab = table(te$Y, pred$Yhat)
  prec_class = diag(result_tab)/colSums(result_tab)
  prec_class[is.na(prec_class)] = 0
  precision =  mean(prec_class)

  recal_class = diag(result_tab)/rowSums(result_tab)
  recall = mean(recal_class)
  f1_score = 2 * (precision * recall) / (precision + recall)

  return(list(precision = precision, recall = recall, f1_score = f1_score))
}
