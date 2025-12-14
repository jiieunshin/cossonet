cossonet.cox = function (x, time, status, nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, cv, nfold, kparam, one.std, scale)
{
  n = length(time)
  p = length(wt)

  message("fit COSSO  with n = ", n, ", p =", ncol(x), "\n")

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
  mscale = rep(1, d)
  message("kernel:", type, "and d =", d, "\n")

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  par(mfrow = c(1,2))
  # solve c (1st)
  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id,  rep(1, d)/mscale^2, lambda0, type, cv, nfold, kparam, one.std = one.std, show = TRUE)

  # solve theta (1st)
  theta_cvfit = cv.gettheta.subset(getc_cvfit, K, time, status, nbasis, basis.id, mscale, getc_cvfit$optlambda, lambda_theta, gamma, cv, nfold, one.std = one.std)

  par(op)

  # solve c (2nd)
  theta.new = rescale_theta(theta_cvfit$theta.new)

  getc_cvfit = cv.getc.subset(K, time, status, nbasis, basis.id, theta.new/mscale^2, lambda0, type, cv, nfold, kparam, one.std = FALSE, show = FALSE)

  out = list(data = list(x = x, time = time, status = status, coord = K$coord, basis.id = basis.id, RS = getc_cvfit$RS, wt = mscale, kernel = type, nfold, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = getc_cvfit,
             theta_step = theta_cvfit,
             family = "Cox")

  return(out)
}

