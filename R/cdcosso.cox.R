#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param x Explanation variable matrix or data frame.
#' @param time Dependent variable vector or matrix or data frame containing time and status columns (for Cox model).
#' @param status Type of statistical model. Use one of the following strings: "gaussian", "binomial", "poisson", "negbin", "svm", or "Cox".
#' @param wt Type of statistical model. Use one of the following strings: "gaussian", "binomial", "poisson", "negbin", "svm", or "Cox".
#' @param lambda0 Type of kernel function to use in case of SVM model. Use one of the following strings: "linear", "gaussian", "poly", "spline", "anova_gaussian", or "gaussian2".
#' @param lambda_theta Type of optimization algorithm. Use either the string "CD" (Coordinate Descent) or "QP".
#' @param gamma Kernel parameter values to use for SVM models.
#' @param nfolds Vector of penalty parameters to be applied to different parts of the model.
#' @param one.std Vector of Lagrange multiplier.
#' @param type Gamma value used in Stochastic Search Optimization.
#' @param kparam Number of folds for cross-validation.
#' @param algo Logical value indicating whether to standardize explanatory variables.
#'
#' @return A list containing information about the fitted model. Depending on the type of dependent variable, various information may be returned.
#' @export

# time = unlist(y[, 'time'])
# stauts = unlist(y[, 'status'])
cdcosso.cox = function (x, time, status, wt, lambda0, lambda_theta, gamma, nfolds, one.std, type, kparam, algo)
{
  # library(survival)
  n = nrow(x)
  d = length(wt)

  par(mfrow = c(3,2))
  # initialize

  # solve theta
  getc_cvfit = cv.getc(x, time, status, rep(1, d)/wt^2, nfolds, lambda0, one.std, type, kparam, algo) ## 초기??? ??????. ????????? ??????
  theta_cvfit = cv.gettheta(getc_cvfit, x, time, status, wt, getc_cvfit$optlambda, lambda_theta, gamma, nfolds, one.std, type, kparam, algo)

  # solve (theta) - 2nd
  theta.new = rescale_theta(theta_cvfit$theta.new, FALSE)
  sspline_cvfit = cv.getc(x, time, status, theta.new/wt^2, nfolds, lambda0, one.std, type, kparam, algo)

  par(mfrow = c(1,1))
  object = list()
  object$family = "Cox"
  out = list(data = list(x = x, time = time, status = status, R = getc_cvfit$R, kernel = type, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = getc_cvfit,
             theta_step = theta_cvfit,
             object,
             algorithm = algo)

  return(out)
}
