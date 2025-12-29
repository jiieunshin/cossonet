cossonet.cox = function (x, time, status, nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, cv, nfold, kparam, one.std, scale)
{
  n = length(time)
  p = length(wt)

  message("fit COSSO  with n = ", n, ", p =", ncol(x), "\n")

  ## --------------------------
  ## (A) scaling step
  ## --------------------------
  x_scaled = apply(x, 2, rescale)
  
  ## --------------------------
  ## (B) basis selection
  ## --------------------------
  if (missing(nbasis) & missing(basis.id)) {
    nbasis = max(40, ceiling(12 * n^(2/9)))
    basis.id = sort(sample(1:n, nbasis))
  }
  if (missing(nbasis) & !missing(basis.id))
    nbasis <- length(basis.id)
  if (!missing(nbasis) & missing(basis.id))
    basis.id <- sort(sample(1:n, nbasis))
  
  nbasis = as.integer(nbasis)
  
  ## --------------------------
  ## (C) kernel construction on TRAIN ONLY
  ## --------------------------
  ## IMPORTANT: use x_scaled (train scale)
  K = make_anovaKernel(x_scaled, x_scaled, type = type, kparam = kparam)
  d = K$numK
  mscale = rep(1, d)
  
  message("kernel: ", type, " and d =", d, "\n")
  
  ## --------------------------
  ## (D) CV for theta (1st step)
  ## --------------------------
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  # par(mfrow = c(1,2))
  
  getc_cvfit = cv.getc.subset(
    K       = K,
    time    = time,
    status  = status,
    nbasis  = nbasis,
    basis.id = basis.id,
    mscale  = rep(1, d)/mscale^2,
    c.init  = NULL,
    cand.lambda = lambda0,
    type    = type,
    cv      = cv,
    nfold   = nfold,
    one.std = FALSE,
    show    = TRUE
  )

  ## --------------------------
  ## (E) CV for (b, c) (2nd step)
  ## --------------------------
  theta_cvfit = cv.gettheta.subset(
    getc_cvfit,
    K        = K,
    time     = time,
    status   = status,
    nbasis   = nbasis,
    basis.id = basis.id,
    mscale   = rep(1, d)/mscale^2,
    lambda0 = getc_cvfit$opt_lambda0,
    lambda_theta = lambda_theta,
    gamma    = gamma,
    cv       = cv,
    nfold    = nfold,
    one.std  = one.std
  )

  theta.new = rescale_theta(theta_cvfit$theta.new)
  
  ## --------------------------
  ## (F) second theta refinement
  ## --------------------------
  getc_cvfit2 = cv.getc.subset(
    K       = K,
    time    = time,
    status  = status,
    nbasis  = nbasis,
    basis.id = basis.id,
    mscale  = theta.new/mscale^2,
    c.init  = getc_cvfit$c.new,
    cand.lambda = getc_cvfit$opt_lambda0,
    type     = type,
    cv       = "mse",
    nfold    = 1,
    one.std  = FALSE,
    show     = FALSE
  )
  
  rm(getc_cvfit)
  # par(mfrow = c(1,1))

  ## --------------------------
  ## (G) save structure
  ## --------------------------
  out = list(
    data = list(
      x         = x_scaled,
      time      = time,
      status    = status,
      coord     = K$coord,
      basis.id  = basis.id,
      ristset   = getc_cvfit$RS,
      wt        = mscale,
      kernel    = type,
      nfold     = nfold,
      kparam    = kparam
    ),
    tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
    c_step = getc_cvfit2,
    theta_step = theta_cvfit,
    family = "Cox"
  )
  
  return(out)
}