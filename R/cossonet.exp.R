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
                         type, cv, nfold, kparam, one.std, scale = TRUE)
{
  n = length(y)
  
  ## --------------------------
  ## (A) scaling step
  ## --------------------------
  if(scale){
    x_mean <- apply(x, 2, mean)
    x_sd   <- apply(x, 2, sd)
    x_sd[x_sd == 0] <- 1
    x_scaled <- scale(x, center = x_mean, scale = x_sd)
  } else {
    x_scaled <- x
    x_mean <- rep(0, ncol(x))
    x_sd   <- rep(1, ncol(x))
  }
  
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
  K = make_anovaKernel(x_scaled, x_scaled[basis.id, ], 
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
    theta0  = rep(1, d)/mscale^2,
    cand.lambda = lambda0,
    obj     = obj,
    type    = type,
    cv      = cv,
    nfold   = nfold,
    kparam  = kparam,
    one.std = one.std,
    show    = TRUE,
    scale.info = list(mean = x_mean, sd = x_sd)   ## pass scaling info
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
    mscale   = mscale,
    optlambda = sspline_cvfit$optlambda,
    cand.lambda = lambda_theta,
    gamma    = gamma,
    cv       = cv,
    nfold    = nfold,
    one.std  = one.std,
    obj      = obj,
    scale.info = list(mean = x_mean, sd = x_sd)
  )
  
  theta.new = rescale_theta(nng_fit$theta.new)
  
  ## --------------------------
  ## (F) second theta refinement
  ## --------------------------
  sspline_cvfit2 = cv.sspline.subset(
    K        = K,
    y        = y,
    nbasis   = nbasis,
    basis.id = basis.id,
    theta0   = theta.new/mscale^2,
    cand.lambda = lambda0,
    obj      = obj,
    type     = type,
    cv       = cv,
    nfold    = nfold,
    kparam   = kparam,
    one.std  = FALSE,
    show     = FALSE,
    scale.info = list(mean = x_mean, sd = x_sd)
  )
  
  ## --------------------------
  ## (G) save structure
  ## --------------------------
  out = list(
    data = list(
      x         = x_scaled,
      raw.x     = x,
      mean      = x_mean,
      sd        = x_sd,
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

