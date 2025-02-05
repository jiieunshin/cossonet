#' The cossonet.cox function is a COSSO model based on the Cox proportional hazards model applied to survival responses.
#' Survival responses should be in matrix form with time and status as column names. It works internally via the cossonet function.
#'
#' @param x Input matrix or data frame of $n$ by $p$. `x` must have at least two columns ($p>1$).
#' @param time for right censored data, this is the follow up time. For interval data, the first argument is the starting time for the interval. It follows the same input format as the `time` argument in the `Surv` function from the `survival` package.
#' @param status Status indicator, typically coded as 0 for "alive" and 1 for "dead." It follows the same input format as the `status` argument in the `Surv` function from the `survival` package.
#' @param nbasis The number of "knots". If `basis.id` is provided, it is set to the length of `basis.id`.
#' @param basis.id The index of the "knot" to select.
#' @param wt The weights assigned to the explanatory variables.
#' @param lambda0 A vector of `lambda0` sequences.
#' @param lambda_theta A vector of `lambda` sequences.
#' @param gamma Elastic-net mixing parameter.
#' @param type The kernel function.
#' @param kparam The kernel function.
#' @param scale Boolean for whether to scale continuous explanatory variables to values between 0 and 1.
#'
#' @return A list of outputs obtained from the fitted model for the Cox proportional hazards model
#' @export

cossonet.cox = function (x, time, status, nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, kparam, scale)
{
  n = length(time)
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
  # solve c (1st)
  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id, rep(1, d)/wt^2, lambda0, type, kparam, one.std = TRUE, show = TRUE)

  # solve theta (1st)
  theta_cvfit = cv.gettheta.subset(getc_cvfit, K, time, status, nbasis, basis.id, wt, getc_cvfit$optlambda, lambda_theta, gamma)

  # solve c (2nd)
  theta.new = rescale_theta(theta_cvfit$theta.new)

  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id, theta.new/wt^2, lambda0, type, kparam, one.std = FALSE, show = FALSE)

  par(mfrow = c(1,1))

  out = list(data = list(x = x, time = time, status = status, basis.id = basis.id, RS = getc_cvfit$RS, wt = wt, kernel = type, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = getc_cvfit,
             theta_step = theta_cvfit,
             family = "Cox")

  return(out)
}

