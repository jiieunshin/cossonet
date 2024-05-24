#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param true A vector of true value of binary class response.
#' @param est A vector of estimated response by `cdcosso`.
#'
#' @return A list of contingency table for predicted results of binary class response.
#' @export
#'

metric = function(true, est){
  result_tab = table(true, est)
  is.col = colnames(result_tab) == c("0", "1")

  # if no one class are exist,
  if(sum(!is.col) > 0){
    colid = which(!is.col)

    if(colid == 1){
      result_tab = cbind(0, result_tab)
      colnames(result_tab) = c("0", "1")
    }

    if(colid == 2){
      result_tab = cbind(result_tab, 0)
      colnames(result_tab) = c("0", "1")
    }
  }

  tp = result_tab[4]
  fp = result_tab[3]
  fn = result_tab[2]
  precision = tp/(tp + fp + 1e-10)
  recall = tp/(tp + fn + 1e-10)

  if((precision + recall) == 0){
    f1_score = 0
  } else{
    f1_score = 2 * (precision * recall)/(precision + recall)
  }

  return(list(tp = tp, fp = fp, precision = precision, recall = recall, f1_score = f1_score))
}

