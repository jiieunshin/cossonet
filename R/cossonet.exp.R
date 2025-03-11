cossonet.exp = function (x, y, wt, nbasis, basis.id, lambda0, lambda_theta, gamma, obj, type, nfold, kparam, one.std, scale)
{
  n = length(y)
  p = length(wt)

  message("fit COSSO  with n = ", n, "p =", ncol(x), "\n")

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
  message("kernel:", type, "and d =", d, "\n")

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  par(mfrow = c(1,2))
  # solve (theta) - 1st
  sspline_cvfit = cv.sspline.subset(K, y, nbasis, basis.id, rep(1, p)/wt^2, lambda0, obj, type, nfold, kparam, one.std = one.std, show = TRUE)

  # solve (b, c) - 1st
  nng_fit = cv.nng.subset(sspline_cvfit, K, y, nbasis, basis.id, wt, sspline_cvfit$optlambda, lambda_theta, gamma, nfold, one.std = one.std, obj)
  theta.new = rescale_theta(nng_fit$theta.new)

  par(op)

  # solve (theta) - 2nd
  sspline_cvfit = cv.sspline.subset(K, y, nbasis, basis.id, rep(1, p) / wt^2, lambda0, obj, type, nfold, kparam, one.std = FALSE, show = FALSE)

  out = list(data = list(x = x, y = y, basis.id = basis.id, wt = wt, kernel = type, kparam = kparam),
             tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
             c_step = sspline_cvfit,
             theta_step = nng_fit,
             family = obj$family)

  return(out)
}

