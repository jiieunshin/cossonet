#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param x Explanation variable matrix or data frame.
#' @param y Dependent variable vector or matrix or data frame containing time and status columns (for Cox model).
#' @param wt Type of statistical model. Use one of the following strings: "gaussian", "binomial", "poisson", "negbin", "svm", or "Cox".
#' @param lambda0 Type of kernel function to use in case of SVM model. Use one of the following strings: "linear", "gaussian", "poly", "spline", "anova_gaussian", or "gaussian2".
#' @param lambda_theta Type of optimization algorithm. Use either the string "CD" (Coordinate Descent) or "QP".
#' @param M Weight vector for each explanatory variable. The default is to use the same weight of 1 for all variables.
#' @param gamma Kernel parameter values to use for SVM models.
#' @param obj Penalty parameter vector for Lasso and Ridge regression.
#' @param nfolds Vector of penalty parameters to be applied to different parts of the model.
#' @param one.std Vector of Lagrange multiplier.
#' @param type Gamma value used in Stochastic Search Optimization.
#' @param kparam Number of folds for cross-validation.
#' @param algo Logical value indicating whether to standardize explanatory variables.
#'
#' @return A list containing information about the fitted model. Depending on the type of dependent variable, various information may be returned.
#' @export

cdcosso.glm = function (x, y, wt, lambda0, lambda_theta, M, gamma, obj, nfolds, one.std, type, kparam, algo)
{
  n = length(y)
  d = length(wt)
  par(mfrow = c(1,2))

  # initiation
  sdx <- sqrt(drop(rep(1, n) %*% (x^2))/(n - 1))
  # init.theta = (as.vector(glmnet(x, y, family = "binomial", lambda = lambda_theta[10], gamma = 0)$beta) + .1) / (sdx + .1)
  init.theta = as.vector(glmnet(x, y, family = "binomial", lambda = lambda_theta[2], gamma = 0)$beta)

  if(sum(init.theta == 0) == d){
    init.theta = rep(1e-10, d)
  } else{
    init.theta = c(scale(init.theta))
  }

  # solve (theta) - 1st
  sspline_cvfit = cv.sspline.cd(x, y, init.theta/wt^2, nfolds, lambda0, obj, one.std, type, kparam, algo) ## 초기값 설정. 수정할 함수
  optlambda0 = sspline_cvfit$optlambda

  # solve (b, c) - 1st
  nng_fit = cv.nng.cd(sspline_cvfit, x, y, wt, init.theta, optlambda0, lambda_theta, M, gamma, nfolds, obj, one.std, algo)

  # solve (theta) - 2nd
  if(sum(nng_fit$theta.new == 0) == d){
    theta.new = rep(1e-10, d)
  } else{
    theta.new = nng_fit$theta.new
  }

  # if(algo == "CD")
  #   Rtheta <- wsGram(sspline_cvfit$R, scale(theta.new)/wt^2)
  #
  # if(algo == "QP")
  #   Rtheta <- wsGram(sspline_cvfit$R, theta.new/wt^2)
  Rtheta <- wsGram(sspline_cvfit$R, theta.new/wt^2)
  sspline_cvfit = cv.sspline.cd(x, y, theta.new/wt^2, nfolds, lambda0, obj, one.std, type, kparam, algo) ## 초기값 설정. 수정할 함수

  par(mfrow = c(1,1))

  if(algo == "CD")
    out = list(data = list(x = x, y = y, R = sspline_cvfit$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
               c_step = sspline_cvfit,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

  if(algo == "QP")
    out = list(data = list(x = x, y = y, R = sspline_cvfit$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, M = M, gamma = gamma),
               c_step = sspline_cvfit,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

  class(out) = "cosso"
  return(out)
}
