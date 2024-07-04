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
#' @param type Gamma value used in Stochastic Search Optimization.
#' @param kparam Number of folds for cross-validation.
#' @param scale Boolean for whether to scale the input data to range between 0 and 1.
#' @param algo Logical value indicating whether to standardize explanatory variables.
#'
#' @return A list containing information about the fitted model. Depending on the type of dependent variable, various information may be returned.
#' @export

# x = tr_x
# time = unlist(tr_y[, "time"])
# status = unlist(tr_y[, "status"])
# type = "spline"
# algo = "CD"
# family = 'Cox'
# gamma = 0.95
# kparam=1
# lambda0 = exp(seq(log(2^{-27}), log(2^{2}), length.out = 20))
# lambda_theta = exp(seq(log(2^{-27}), log(2^{2}), length.out = 20))
# wt = rep(1, ncol(x))
cdcosso.cox = function (x, time, status, wt, lambda0, lambda_theta, gamma, type, kparam, scale, algo)
{
  n = length(time)
  p = length(wt)
  # cat("fit COSSO  with n = ", n, "p =", p, "\n")

  K = make_anovaKernel(x, x, type = type, kparam, scale)
  d = K$numK
  # cat("kernel:", type, "and d =", d, "\n")

  par(mfrow = c(1,3))
  # solve theta
  getc_cvfit  = cv.getc(K, time, status, rep(1, d)/wt^2, lambda0, type, kparam, algo, show = TRUE)
  theta_cvfit = cv.gettheta(getc_cvfit, x, time, status, wt, getc_cvfit$optlambda, lambda_theta, gamma, type, kparam, algo)

  # solve (theta) - 2nd
  theta.new = rescale_theta(theta_cvfit$theta.new)
  # print(theta.new)
  getc_cvfit = cv.getc(K, time, status, theta.new/wt^2, lambda0, type, kparam, algo, show = TRUE)

  par(mfrow = c(1,1))

  out = list(data = list(x = x, time = time, status = status, RiskSet = RiskSet(time, status), R = getc_cvfit$R, kernel = type, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = getc_cvfit,
             theta_step = theta_cvfit,
             family = "Cox",
             algorithm = algo)

  return(out)
}
