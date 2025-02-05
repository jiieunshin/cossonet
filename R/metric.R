#' The function `metric` provides a contingency table for the predicted class and the true class for binary classes.
#'
#' @param true binary true class.
#' @param est binary predicted class.
#'
#' @return a contingency table for the predicted results of binary class responses.
#' @export
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

