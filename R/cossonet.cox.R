cossonet.cox = function (x, time, status, nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, nfold, kparam, one.std, scale)
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
  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id, rep(1, d)/wt^2, lambda0, type, nfold, kparam, one.std = one.std, show = TRUE)

  # solve theta (1st)
  theta_cvfit = cv.gettheta.subset(getc_cvfit, K, time, status, nbasis, basis.id, wt, getc_cvfit$optlambda, lambda_theta, gamma, nfold, one.std = one.std)

  # solve c (2nd)
  theta.new = rescale_theta(theta_cvfit$theta.new)

  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id, theta.new/wt^2, lambda0, type, nfold, kparam, one.std = FALSE, show = FALSE)

  par(mfrow = c(1,1))

  out = list(data = list(x = x, time = time, status = status, basis.id = basis.id, RS = getc_cvfit$RS, wt = wt, kernel = type, nfold, kparam = kparam, one.std = one.std),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = getc_cvfit,
             theta_step = theta_cvfit,
             family = "Cox")

  return(out)
}

