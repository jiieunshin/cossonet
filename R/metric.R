#' The function `metric` provides a contingency table for the predicted class and the true class for binary classes.
#'
#' @param true binary true class.
#' @param est binary predicted class.
#'
#' @return a contingency table for the predicted results of binary class responses.
#'
#' @examples
#' \donttest{
#' set.seed(20250101)
#' tr = data_generation(n = 200, p = 20, SNR = 9, response = "continuous")
#' tr_x = tr$x
#' tr_y = tr$y
#'
#' te = data_generation(n = 1000, p = 20, SNR = 9, response = "continuous")
#' te_x = te$x
#' te_y = te$y
#'
#' # Fit the model
#' fit = cossonet(tr_x, tr_y, family = 'gaussian', gamma = 0.95, kernel = "spline", scale = TRUE,
#'       lambda0 = exp(seq(log(2^{-4}), log(2^{0}), length.out = 20)),
#'       lambda_theta = exp(seq(log(2^{-8}), log(2^{-6}), length.out = 20))
#'       )
#'
#' # Predict new dataset
#' pred = cossonet.predict(fit, te_x)
#'
#' # Calculate the contingency table for binary class
#' true_var = c(rep(1, 4), rep(0, 20-4))
#' est_var = ifelse(fit$theta_step$theta.new > 0, 1, 0)
#' metric(true_var, est_var)
#' }
#'
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

