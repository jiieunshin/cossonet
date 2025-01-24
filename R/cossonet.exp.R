#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param x Input matrix with $n$ observations and $p$ dimension.
#' @param y Response variable. The type can be continuous (default), binary class, non-negative count, survival.
#' @param wt A weight vector for components.
#' @param nbasis The number of basis.
#' @param basis.id The index of basis.
#' @param lambda0 A vector of tuning parameter to select optimal smoothing parameter.
#' @param lambda_theta A vector of tuning parameter to select optimal tuning parameter.
#' @param gamma The elastic net mixing parameter between 0 and 1. When `gamma = 1`, the penalty is to be LASSO. When `gamma = 0`, the penalty is to be ridge penalty. The default is `gamma = 0.95`.
#' @param obj The type of family.
#' @param type Kernel function which is used to convert the input data for training and predicting. The four types is provided, `linear` (default), `gaussian`, `poly`, and `spline`.
#' @param kparam Kernel parameter values that is used in gaussian kernel and polynomial kernel.
#' @param scale Boolean for whether to scale the input data to range between 0 and 1.
#'
#' @return A list containing information about the fitted model. Depending on the type of dependent variable, various information may be returned.
#' @export

cossonet.exp = function (x, y, wt, nbasis, basis.id, lambda0, lambda_theta, gamma, obj, type, kparam, scale)
{
  n = length(y)
  p = length(wt)

  cat("fit COSSO  with n = ", n, "p =", ncol(x), "\n")

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
  cat("kernel:", type, "and d =", d, "\n")

  par(mfrow = c(1,2))
  # solve (theta) - 1st
  sspline_cvfit = cv.sspline.subset(K, y, nbasis, basis.id, rep(1, p)/wt^2, lambda0, obj, type, kparam, one.std = TRUE,  show = TRUE)

  # solve (b, c) - 1st
  nng_fit = cv.nng.subset(sspline_cvfit, K, y, nbasis, basis.id, wt, sspline_cvfit$optlambda, lambda_theta, gamma, obj)
  theta.new = rescale_theta(nng_fit$theta.new)

  # solve (theta) - 2nd
  sspline_cvfit = try({cv.sspline.subset(K, y, nbasis, basis.id, rep(1, p)/wt^2, lambda0, obj, type, kparam, one.std = FALSE, show = FALSE)})

  par(mfrow = c(1,1))

  out = list(data = list(x = x, y = y, Uv = sspline_cvfit$Uv, basis.id = basis.id, wt = wt, kernel = type, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = sspline_cvfit,
             theta_step = nng_fit,
             family = obj$family)

  return(out)
}

