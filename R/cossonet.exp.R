# cossonet.exp = function (x, y, wt, nbasis, basis.id, lambda0, lambda_theta, gamma, obj, type, cv, nfold, kparam, one.std, scale)
# {
#   n = length(y)
#   p = length(wt)
# 
#   message("fit COSSO  with n = ", n, ", p =", ncol(x), "\n")
# 
#   if (missing(nbasis) & missing(basis.id)) {
#     nbasis = max(40, ceiling(12 * n^(2/9)))
#     basis.id = sort(sample(1:n, nbasis))
#   }
#   if (missing(nbasis) & !missing(basis.id))
#     nbasis <- length(basis.id)
#   if (!missing(nbasis) & missing(basis.id))
#     basis.id <- sort(sample(1:n, nbasis))
# 
#   nbasis = as.integer(nbasis)
# 
#   K = make_anovaKernel(x, x, type = type, kparam, scale)
#   d = K$numK
#   mscale = rep(1, d)
# 
#   message("kernel: ", type, " and d =", d, "\n")
# 
#   op <- par(no.readonly = TRUE)
#   on.exit(par(op))
# 
#   par(mfrow = c(1,2))
#   # solve (theta) - 1st
#   sspline_cvfit = cv.sspline.subset(K, y, nbasis, basis.id, rep(1, d)/mscale^2, lambda0, obj, type, cv, nfold, kparam, one.std = one.std, show = TRUE)
# 
#   # solve (b, c) - 1st
#   nng_fit = cv.nng.subset(sspline_cvfit, K, y, nbasis, basis.id, mscale, sspline_cvfit$optlambda, lambda_theta, gamma, cv, nfold, one.std = one.std, obj)
#   theta.new = rescale_theta(nng_fit$theta.new)
# 
#   # solve (theta) - 2nd
#   sspline_cvfit = cv.sspline.subset(K, y, nbasis, basis.id, theta.new/mscale^2, lambda0, obj, type, cv, nfold, kparam, one.std = FALSE, show = FALSE)
#   # nng_fit = cv.nng.subset(sspline_cvfit, K, y, nbasis, basis.id, mscale, sspline_cvfit$optlambda, lambda_theta, gamma, cv, nfold, one.std = one.std, obj)
#   
#   out = list(data = list(x = x, y = y, coord = K$coord, basis.id = basis.id, wt = mscale, kernel = type, kparam = kparam),
#              tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
#              c_step = sspline_cvfit,
#              theta_step = nng_fit,
#              family = obj$family)
# 
#   return(out)
# }

cossonet.exp = function (x, y, wt, nbasis, basis.id, 
                         lambda0, lambda_theta, gamma, obj, 
                         type, cv, nfold, kparam, one.std, scale)
{
  n = length(y)
  
  ## --------------------------
  ## (A) scaling step
  ## --------------------------
  x_scaled = apply(x, 2, rescale)
  
  ## weight dimension check
  p = length(wt)
  
  message("fit COSSO  with n = ", n, ", p =", ncol(x), "\n")
  
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
  K = make_anovaKernel(x_scaled, x_scaled, 
                       type = type, kparam = kparam)
  d = K$numK
  mscale = rep(1, d)
  
  message("kernel: ", type, " and d =", d, "\n")
  
  ## --------------------------
  ## (D) CV for theta (1st step)
  ## --------------------------
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(1,2))
  
  sspline_cvfit = cv.sspline.subset(
    K       = K,
    y       = y,
    nbasis  = nbasis,
    basis.id = basis.id,
    mscale  = rep(1, d)/mscale^2,
    c.init  = NULL,
    cand.lambda = lambda0,
    obj     = obj,
    type    = type,
    cv      = "GCV",
    nfold   = nfold,
    one.std = one.std,
    show    = TRUE
  )
  
  ## --------------------------
  ## (E) CV for (b, c) (2nd step)
  ## --------------------------
  nng_fit = cv.nng.subset(
    sspline_cvfit,
    K        = K,
    y        = y,
    nbasis   = nbasis,
    basis.id = basis.id,
    mscale   = rep(1, d)/mscale^2,
    lambda0 = sspline_cvfit$opt_lambda0,
    lambda_theta = lambda_theta,
    gamma    = gamma,
    cv       = cv,
    nfold    = nfold,
    one.std  = one.std,
    obj      = obj
  )
  
  # theta.new = rescale_theta(nng_fit$theta.new)
  
  ## --------------------------
  ## (F) second theta refinement
  ## --------------------------
  sspline_cvfit2 = cv.sspline.subset(
    K        = K,
    y        = y,
    nbasis   = nbasis,
    basis.id = basis.id,
    mscale   = nng_fit$theta.new/mscale^2,
    c.init   = sspline_cvfit$c.new,
    cand.lambda = sspline_cvfit$opt_lambda0,
    obj      = obj,
    type     = type,
    cv       = "GCV",
    nfold    = 1,
    one.std  = one.std,
    show     = FALSE
  )
  
  rm(sspline_cvfit)
  
  ## --------------------------
  ## (G) save structure
  ## --------------------------
  out = list(
    data = list(
      x         = x_scaled,
      y         = y,
      coord     = K$coord,
      basis.id  = basis.id,
      wt        = mscale,
      kernel    = type,
      kparam    = kparam
    ),
    tune = list(lambda0 = lambda0, lambda_theta = lambda_theta, gamma = gamma),
    c_step = sspline_cvfit2,
    theta_step = nng_fit,
    family = obj$family
  )
  
  return(out)
}

