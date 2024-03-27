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

cdcosso.glm = function (x, y, wt, lambda0, lambda_theta, gamma, obj, nfolds, one.std, type, kparam, algo)
{
  n = length(y)
  d = length(wt)
  par(mfrow = c(3,2))

  # initiation
  # init.theta = as.vector(glmnet(x, y, family = "binomial", lambda = lambda_theta[2], gamma = 0)$beta)

  # solve (theta) - 1st
  sspline_cvfit1 = cv.sspline(x, y, rep(1, d)/wt^2, nfolds, lambda0, obj, one.std, type, kparam, algo) ## 초기값 설정. 수정할 함수

  # solve (b, c) - 1st
  nng_fit = cv.nng(sspline_cvfit1, x, y, wt, sspline_cvfit1$optlambda, lambda_theta, gamma, nfolds, obj, one.std, algo)
  theta.new = rescale_theta(nng_fit$theta.new, FALSE)

  # solve (theta) - 2nd
  sspline_cvfit2 = try({cv.sspline(x, y, theta.new/wt^2, nfolds, lambda0, obj, one.std, type, kparam, algo)}) ## 초기값 설정. 수정할 함수

  # if not convergence
  if(class(sspline_cvfit2) == "try-error"){
    out = list(data = list(x = x, y = y, R = sspline_cvfit1$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
               c_step = sspline_cvfit1,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

    class(out) = "cosso"
    cat("cdcosso is not convergence. \n")
    return(out)
  } else if(!sspline_cvfit2$conv){
    out = list(data = list(x = x, y = y, R = sspline_cvfit1$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
               c_step = sspline_cvfit1,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

    class(out) = "cosso"
    cat("cdcosso is not convergence. \n")
    return(out)
  }

  nng_fit = cv.nng(sspline_cvfit2, x, y, wt, sspline_cvfit2$optlambda, lambda_theta, gamma, nfolds, obj, one.std, algo)
  theta.new = rescale_theta(nng_fit$theta.new, FALSE)

  print(nng_fit$theta.new)
  sspline_cvfit3 = cv.sspline(x, y, theta.new/wt^2, nfolds, lambda0, obj, one.std, type, kparam, algo) ## 초기값 설정. 수정할 함수


  # if not convergence
  if(class(sspline_cvfit3) == "try-error"){
    out = list(data = list(x = x, y = y, R = sspline_cvfit2$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
               c_step = sspline_cvfit2,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

    class(out) = "cosso"
    cat("cdcosso is not convergence. \n")
    return(out)
  } else if(!sspline_cvfit3$conv){
    out = list(data = list(x = x, y = y, R = sspline_cvfit2$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
               c_step = sspline_cvfit2,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

    class(out) = "cosso"
    cat("cdcosso is not convergence. \n")
    return(out)
  }


  if(algo == "CD")
    out = list(data = list(x = x, y = y, R = sspline_cvfit3$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
               c_step = sspline_cvfit3,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

  if(algo == "QP")
    out = list(data = list(x = x, y = y, R = sspline_cvfit3$R, kernel = type, kparam = kparam),
               tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
               c_step = sspline_cvfit3,
               theta_step = nng_fit,
               object = obj,
               algorithm = algo)

  class(out) = "cosso"
  return(out)
}
