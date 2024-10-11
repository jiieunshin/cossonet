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
#' @param nbasis The number of basis.
#' @param basis.id The index of basis.
#' @param wt Type of statistical model. Use one of the following strings: "gaussian", "binomial", "poisson", "negbin", "svm", or "Cox".
#' @param lambda0 Type of kernel function to use in case of SVM model. Use one of the following strings: "linear", "gaussian", "poly", "spline", "anova_gaussian", or "gaussian2".
#' @param lambda_theta Type of optimization algorithm. Use either the string "CD" (Coordinate Descent) or "QP".
#' @param gamma Kernel parameter values to use for SVM models.
#' @param type Gamma value used in Stochastic Search Optimization.
#' @param kparam Number of folds for cross-validation.
#' @param scale Boolean for whether to scale the input data to range between 0 and 1.
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
# lambda0 = exp(seq(log(2^{-8}), log(2^{-6}), length.out = 20))
# lambda_theta = exp(seq(log(2^{-8}), log(2^{-6}), length.out = 20))
# wt = rep(1, ncol(x))
cdcosso.cox = function (x, time, status, nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, kparam, scale)
{
  n = length(time)
  p = length(wt)
  # cat("fit COSSO  with n = ", n, "p =", p, "\n")

  if (missing(nbasis) & missing(basis.id)) {
    nbasis = max(40, ceiling(12 * n^(2/9)))
    basis.id = sort(sample(1:n, nbasis))
  }
  if (missing(nbasis) & !missing(basis.id))
    nbasis <- length(basis.id)
  if (!missing(nbasis) & missing(basis.id))
    basis.id <- sort(sample(1:n, nbasis))

  nbasis = as.integer(nbasis)

  K = make_anovaKernel(x, x, type = type, kparam, scale)
  d = K$numK
  # cat("kernel:", type, "and d =", d, "\n")

  par(mfrow = c(1,3))
  # solve c (1st)
  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id, rep(1, d)/wt^2, lambda0, type, kparam, one.std = TRUE, show = TRUE)

  # solve theta (1st)
  theta_cvfit = cv.gettheta.subset(getc_cvfit, K, time, status, nbasis, basis.id, wt, getc_cvfit$optlambda, lambda_theta, gamma)

  # solve c (2nd)
  theta.new = rescale_theta(theta_cvfit$theta.new)
  # print(theta.new)
  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id, theta.new/wt^2, lambda0, type, kparam, one.std = FALSE, show = TRUE)

  # solve theta (2nd)
  # theta_cvfit = cv.gettheta(getc_cvfit, x, time, status, wt, getc_cvfit$optlambda, lambda_theta, gamma, type, kparam)

  par(mfrow = c(1,1))

  out = list(data = list(x = x, time = time, status = status, basis.id = basis.id, R = getc_cvfit$R, RS = getc_cvfit$RS, wt = wt, kernel = type, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = getc_cvfit,
             theta_step = theta_cvfit,
             family = "Cox")

  return(out)
}

