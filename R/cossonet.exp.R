#' The cossonet.exp function is a COSSO model applied to responses belonging to the exponential family.
#' This function is applied to continuous, binary, and count responses.
#' It operates internally through the cossonet function.
#'
#' @param x Input matrix or data frame of $n$ by $p$. `x` must have at least two columns ($p>1$).
#' @param y A response vector with a continuous, binary, or count type.
#' @param wt The weights assigned to the explanatory variables.
#' @param nbasis The number of "knots". If `basis.id` is provided, it is set to the length of `basis.id`.
#' @param basis.id The index of the "knot" to select.
#' @param lambda0 A vector of `lambda0` sequences.
#' @param lambda_theta A vector of `lambda` sequences.
#' @param gamma Elastic-net mixing parameter.
#' @param obj distribution used in the model.
#' @param type The kernel function. One of four types of `linear` (default), `gaussian`, `poly`, and `spline`.
#' @param nfold The number of folds to use in cross-validation is used to determine how many subsets to divide the data into for the training and validation sets.
#' @param kparam Parameters for Gaussian and polynomial kernel functions.
#' @param one.std A logical value indicating whether to apply the "1-standard error rule," selecting the simplest model within one standard error of the best model.
#' @param scale Boolean for whether to scale continuous explanatory variables to values between 0 and 1. Default of `TRUE`.
#'
#' @return A list of outputs obtained from a model fitted to an exponential distribution.
#' @export

cossonet.exp = function (x, y, wt, nbasis, basis.id, lambda0, lambda_theta, gamma, obj, type, nfold, kparam, one.std, scale)
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
  sspline_cvfit = cv.sspline.subset(K, y, nbasis, basis.id, rep(1, p)/wt^2, lambda0, obj, type, nfold, kparam, one.std = one.std, show = TRUE)

  # solve (b, c) - 1st
  nng_fit = cv.nng.subset(sspline_cvfit, K, y, nbasis, basis.id, wt, sspline_cvfit$optlambda, lambda_theta, gamma, nfold, one.std = one.std, obj)
  theta.new = rescale_theta(nng_fit$theta.new)

  # solve (theta) - 2nd
  sspline_cvfit = try({cv.sspline.subset(K, y, nbasis, basis.id, rep(1, p)/wt^2, lambda0, obj, type, nfold, kparam, one.std = FALSE, show = FALSE)})

  par(mfrow = c(1,1))

  out = list(data = list(x = x, y = y, basis.id = basis.id, wt = wt, kernel = type, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = sspline_cvfit,
             theta_step = nng_fit,
             family = obj$family)

  return(out)
}

