# cv.sspline = function (K, y, mscale, cand.lambda, obj, type, kparam, algo, show)
# {
#   cat("-- c-step -- \n")
#   cat("proceeding... \n")
#   d = K$numK
#   n <- length(y)
#   len = length(cand.lambda)
#
#   R = array(NA, c(n, n, d))
#   for(j in 1:d){
#     R[, , j] = K$K[[j]]
#   }
#
#   Rtheta <- combine_kernel(R, mscale)
#   f.init = rep(0.5, n)
#
#   measure <- rep(NA, len)
#   for (k in 1:len) {
#
#     if(algo == "CD"){
#
#       # initialize
#       ff = f.init
#       mu = obj$linkinv(ff)
#       w = obj$variance(mu)
#       z = ff + (y - mu) / w
#
#       c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = cand.lambda[k])$beta)
#       zw = z * sqrt(w)
#       Rw = Rtheta * w
#       cw = c.init / sqrt(w)
#       sw = sqrt(w)
#
#       fit = .Call("glm_c_step", zw, Rw, Rw, cw, sw, n, n, cand.lambda[k], PACKAGE = "cdcosso")
#       # fit = sspline.cd(Rtheta, y, ff, cand.lambda[k], obj, c.init)
#       b.new = fit$b.new
#       cw.new = fit$cw.new
#       c.new = cw.new * sqrt(w)
#     }
#
#     if(sum(is.nan(cw.new)) == n){
#       next
#     } else{
#       ff = f.init
#       mu = obj$linkinv(ff)
#       w = obj$variance(mu)
#       Rw = Rtheta * w
#
#       # validation
#       testfhat = c(b.new + Rtheta %*% c.new)
#       testmu = obj$linkinv(testfhat)
#       testw = obj$variance(testmu)
#
#       XX = zw - Rw %*% cw.new - b.new * sqrt(w)
#       num = t(XX) %*% XX + 1
#       # den = (1 - sum(diag(Rtheta %*% ginv(Rtheta + diag(w)/cand.lambda[k]))) / n)^2
#       S = Rw %*% ginv(t(Rw) %*% Rw) %*% t(Rw)
#       den = (1 - sum(diag(S)) / n)^2 + 1
#       # measure[k] <- as.vector( num / den / n)
#       measure[k] = KL(testfhat, obj)
#     }
#   }
#
#   # plotting error bar
#   if(obj$family == 'gaussian'){
#     main = "Gaussian Family"
#   }
#   if(obj$family == 'binomial'){
#     main = "Binomial Family"
#   }
#   if(obj$family == 'poisson'){
#     main = "Poisson Family"
#   }
#
#   ylab = expression("GCV(" * lambda[0] * ")")
#
#   # optimal lambda1
#   id = which.min(measure)
#   optlambda = cand.lambda[id]
#
#   if(show) plot(log(cand.lambda), measure, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab, ylim = range(measure), pch = 15, col = 'red')
#
#
#   mu = obj$linkinv(f.init)
#   w = obj$variance(mu)
#   z = f.init + (y - mu) / w
#
#   c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = optlambda)$beta)
#   zw = z * sqrt(w)
#   Rw = Rtheta * w
#   cw = c.init / sqrt(w)
#   sw = sqrt(w)
#
#   fit = .Call("glm_c_step", zw, Rw, Rw, cw, sw, n, n, optlambda, PACKAGE = "cdcosso")
#   # fit = sspline.cd(Rtheta, y, ff, optlambda, obj, c.init)
#
#   b.new = fit$b.new
#   cw.new = fit$cw.new
#   c.new = cw.new * sqrt(w)
#
#   f.new = c(b.new + Rtheta %*% c.new)
#   mu.new = obj$linkinv(f.new)
#   w.new = obj$variance(mu.new)
#   z.new = f.new + (y - mu.new) / w.new
#
#   if(obj$family == "binomial") miss <- mean(y != ifelse(mu.new < 0.5, 0, 1))
#   if(obj$family == "gaussian") miss <- mean((y - f.new)^2)
#   if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))
#
#   cat("training error:", miss, "\n")
#
#   out = list(measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
#              z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = b.new,
#              cw.new = cw.new, c.new = c.new, optlambda = optlambda, conv = TRUE)
#
#   rm(K)
#   rm(Rtheta)
#
#   return(out)
# }

# cand.lambda = lambda0
# mscale = wt
# obj = binomial()
cv.sspline.subset = function (K, y, nbasis, basis.id, mscale, cand.lambda, obj, type, kparam, algo, show)
{
  cat("-- c-step -- \n")
  cat("proceeding... \n")
  d = K$numK
  n <- length(y)
  len = length(cand.lambda)

  R = array(NA, c(n, nbasis, d))
  for(j in 1:d){
    R[, , j] = K$K[[j]][, basis.id]
  }

  Rtheta <- combine_kernel(R, mscale)

  R2 = array(NA, c(nbasis, nbasis, d))
  for(j in 1:d){
    R2[, , j] = K$K[[j]][basis.id, basis.id]
  }

  Rtheta2 <- combine_kernel(R2, mscale)

  # initialize
  f.init = rep(0.5, n)
  ff = f.init
  mu = obj$linkinv(ff)
  w = obj$variance(mu)
  z = ff + (y - mu) / w

  #
  fold = cvsplitID(n, 5, y, family = obj$family)
  measure <- matrix(NA, 5, len)
  for(f in 1:5){
    tr_id = as.vector(fold[, -f])
    te_id = fold[, f]

    tr_id = tr_id[!is.na(tr_id)]
    te_id = te_id[!is.na(te_id)]

    tr_n = length(tr_id)
    te_n = length(te_id)

    tr_R = array(NA, c(tr_n, nbasis, d))
    for(j in 1:d){
      tr_R[, , j] = K$K[[j]][tr_id, basis.id]
    }

    tr_Rtheta <- combine_kernel(tr_R, mscale)

    te_R = array(NA, c(te_n, nbasis, d))
    for(j in 1:d){
      te_R[, , j] = K$K[[j]][te_id, basis.id]
    }

    te_Rtheta <- combine_kernel(te_R, mscale)

    #
    zw = z[tr_id] * sqrt(w)[tr_id]
    Rw = tr_Rtheta * w[tr_id]
    sw = sqrt(w)[tr_id]

    for (k in 1:len){

      if(algo == "CD"){
        c.init = as.vector(glmnet(Rw, zw, family = 'gaussian', lambda = cand.lambda[k])$beta)

        # cw = c.init / sqrt(w)[basis.id]

        fit = .Call("glm_c_step", zw, Rw, Rtheta2, c.init, sw, tr_n, nbasis, tr_n * cand.lambda[k], PACKAGE = "cdcosso")
        b.new = fit$b.new
        c.new = fit$cw.new
        # c.new = cw.new * sqrt(w)[basis.id]
        # cat("R calculate:", sum(zw - Rw %*% cw.new) / sum(sw), "\n")
        # cat("C calculate:", b.new, "\n")
      }

      # validation
      testfhat = c(b.new + te_Rtheta %*% c.new)
      testmu = obj$linkinv(testfhat)

      # if(obj$family == "gaussian") measure[f, k] <- mean((testfhat - y[te_id])^2)
      # if(obj$family == "binomial") measure[f, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
      # if(obj$family == "poisson") measure[f, k] <- mean(poisson()$dev.resids(y[te_id], testmu, rep(1, te_n)))
      measure[f, k] <- KL(testfhat, obj)

    }
  }

  # plotting error bar
  if(obj$family == 'gaussian'){
    main = "Gaussian Family"
  }
  if(obj$family == 'binomial'){
    main = "Binomial Family"
  }
  if(obj$family == 'poisson'){
    main = "Poisson Family"
  }

  ylab = expression("GCV(" * lambda[0] * ")")

  # optimal lambda1
  measure_mean = colMeans(measure, na.rm = T)
  measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(te_n)

  sel_id = which(!is.nan(measure_se) & measure_se != Inf)
  measure_mean = measure_mean[sel_id]
  measure_se = measure_se[sel_id]
  cand.lambda = cand.lambda[sel_id]

  min_id = which.min(measure_mean)
  optlambda = cand.lambda[min_id]

  if(show){
    plot(log(cand.lambda), measure_mean, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
         ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
    arrows(x0 = log(cand.lambda), y0 = measure_mean - measure_se,
           x1 = log(cand.lambda), y1 = measure_mean + measure_se,
           angle = 90, code = 3, length = 0.1, col = "darkgray")
  }

  rm(tr_R)
  rm(te_R)
  rm(tr_Rtheta)
  rm(te_Rtheta)

  zw = z * sqrt(w)
  Rw = Rtheta * w
  sw = sqrt(w)

  c.init = as.vector(glmnet(Rw, zw, family = 'gaussian', lambda = optlambda)$beta)

  # cw = c.init / sqrt(w)[basis.id]

  fit = .Call("glm_c_step", zw, Rw, Rtheta2, c.init, sw, n, nbasis, n * optlambda, PACKAGE = "cdcosso")
  b.new = fit$b.new
  c.new = fit$cw.new

  f.new = c(b.new + Rtheta %*% c.new)
  mu.new = obj$linkinv(f.new)
  w.new = obj$variance(mu.new)
  z.new = f.new + (y - mu.new) / w.new

  if(obj$family == "binomial") miss <- mean(y != ifelse(mu.new < 0.5, 0, 1))
  if(obj$family == "gaussian") miss <- mean((y - f.new)^2)
  if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))

  cat("training error:", miss, "\n")

  rm(K)
  rm(Rtheta)
  rm(Rw)

  out = list(measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
             z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = b.new,
             c.new = c.new, optlambda = optlambda, conv = TRUE)


  return(out)
}


# cv.sspline = function (K, y, mscale, cand.lambda, obj, type, kparam, algo, show)
# {
#   cat("-- c-step -- \n")
#   cat("proceeding... \n")
#   d = K$numK
#   n <- length(y)
#   len = length(cand.lambda)
#
#   R = array(NA, c(n, nbasis, d))
#   for(j in 1:d){
#     R[, , j] = K$K[[j]][, basis.id]
#   }
#
#   Rtheta <- combine_kernel(R, mscale)
#
#   R2 = array(NA, c(nbasis, nbasis, d))
#   for(j in 1:d){
#     R2[, , j] = K$K[[j]][basis.id, basis.id]
#   }
#
#   Rtheta2 <- combine_kernel(R2, mscale)
#
#   # initialize
#   f.init = rep(0.5, n)
#   ff = f.init
#   mu = obj$linkinv(ff)
#   w = obj$variance(mu)
#   z = ff + (y - mu) / w
#
#   #
#   fold = cvsplitID(n, 5, y, family = obj$family)
#   measure <- matrix(NA, 5, len)
#   for(f in 1:5){
#     tr_id = as.vector(fold[, -f])
#     te_id = fold[, f]
#
#     tr_id = tr_id[!is.na(tr_id)]
#     te_id = te_id[!is.na(te_id)]
#
#     tr_n = length(tr_id)
#     te_n = length(te_id)
#
#     tr_R = array(NA, c(tr_n, nbasis, d))
#     for(j in 1:d){
#       tr_R[, , j] = K$K[[j]][tr_id, basis.id]
#     }
#
#     tr_Rtheta <- combine_kernel(tr_R, mscale)
#
#     te_R = array(NA, c(te_n, nbasis, d))
#     for(j in 1:d){
#       te_R[, , j] = K$K[[j]][te_id, basis.id]
#     }
#
#     te_Rtheta <- combine_kernel(te_R, mscale)
#
#     #
#     zw = z[tr_id] * sqrt(w)[tr_id]
#     Rw = tr_Rtheta * w[tr_id]
#     sw = sqrt(w)[tr_id]
#
#     for (k in 1:len){
#
#       if(algo == "CD"){
#         c.init = as.vector(glmnet(Rw, zw, family = 'gaussian', lambda = cand.lambda[k])$beta)
#
#         # cw = c.init / sqrt(w)[basis.id]
#
#         fit = .Call("glm_c_step", zw, Rw, Rtheta2, c.init, sw, tr_n, nbasis, tr_n * cand.lambda[k], PACKAGE = "cdcosso")
#         b.new = fit$b.new
#         c.new = fit$cw.new
#         # c.new = cw.new * sqrt(w)[basis.id]
#         # cat("R calculate:", sum(zw - Rw %*% cw.new) / sum(sw), "\n")
#         # cat("C calculate:", b.new, "\n")
#       }
#
#       # validation
#       testfhat = c(b.new + te_Rtheta %*% c.new)
#       testmu = obj$linkinv(testfhat)
#
#       # if(obj$family == "gaussian") measure[f, k] <- mean((testfhat - y[te_id])^2)
#       # if(obj$family == "binomial") measure[f, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
#       # if(obj$family == "poisson") measure[f, k] <- mean(poisson()$dev.resids(y[te_id], testmu, rep(1, te_n)))
#       measure[f, k] <- KL(testfhat, obj)
#
#     }
#   }
#
#   # plotting error bar
#   if(obj$family == 'gaussian'){
#     main = "Gaussian Family"
#   }
#   if(obj$family == 'binomial'){
#     main = "Binomial Family"
#   }
#   if(obj$family == 'poisson'){
#     main = "Poisson Family"
#   }
#
#   ylab = expression("GCV(" * lambda[0] * ")")
#
#   # optimal lambda1
#   measure_mean = colMeans(measure, na.rm = T)
#   measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(te_n)
#
#   sel_id = which(!is.nan(measure_se) & measure_se != Inf)
#   measure_mean = measure_mean[sel_id]
#   measure_se = measure_se[sel_id]
#   cand.lambda = cand.lambda[sel_id]
#
#   min_id = which.min(measure_mean)
#   optlambda = cand.lambda[min_id]
#
#   if(show){
#     plot(log(cand.lambda), measure_mean, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
#          ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
#     arrows(x0 = log(cand.lambda), y0 = measure_mean - measure_se,
#            x1 = log(cand.lambda), y1 = measure_mean + measure_se,
#            angle = 90, code = 3, length = 0.1, col = "darkgray")
#   }
#
#   rm(tr_R)
#   rm(te_R)
#   rm(tr_Rtheta)
#   rm(te_Rtheta)
#
#   zw = z * sqrt(w)
#   Rw = Rtheta * w
#   sw = sqrt(w)
#
#   c.init = as.vector(glmnet(Rw, zw, family = 'gaussian', lambda = optlambda)$beta)
#
#   # cw = c.init / sqrt(w)[basis.id]
#
#   fit = .Call("glm_c_step", zw, Rw, Rtheta2, c.init, sw, n, nbasis, n * optlambda, PACKAGE = "cdcosso")
#   b.new = fit$b.new
#   c.new = fit$cw.new
#
#   f.new = c(b.new + Rtheta %*% c.new)
#   mu.new = obj$linkinv(f.new)
#   w.new = obj$variance(mu.new)
#   z.new = f.new + (y - mu.new) / w.new
#
#   if(obj$family == "binomial") miss <- mean(y != ifelse(mu.new < 0.5, 0, 1))
#   if(obj$family == "gaussian") miss <- mean((y - f.new)^2)
#   if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))
#
#   cat("training error:", miss, "\n")
#
#   rm(K)
#   rm(Rtheta)
#   rm(Rw)
#
#   out = list(measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
#              z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = b.new,
#              c.new = c.new, optlambda = optlambda, conv = TRUE)
#
#
#   return(out)
# }


sspline.cd = function (R, y, f, lambda0, obj, c.init)
{
  n = length(y)
  mu = obj$linkinv(f)

  # initialize
  w = obj$variance(mu)
  z = f + (y - mu) / w
  b = 0

  zw = z * sqrt(w)
  Rw = R * w
  cw = c.init / sqrt(w)
  sw = sqrt(w)
  cw.new = rep(0, n)
  for(i in 1:15){ # outer iteration

    for(j in 1:n){
      L = 2 * sum((zw - Rw[,-j] %*% cw[-j] - b * sw) * Rw[,j]) - n * lambda0 * c(Rw[j,-j] %*% cw[-j])
      R = 2 * sum(Rw[,j]^2) + n * lambda0 * Rw[j,j]
      cw.new[j] = L/R

      loss = abs(cw-cw.new)
      conv1 = max(loss) < 1e-6
      conv2 = min(loss) > 10

      if(conv1 | conv2) break
      cw[j] = cw.new[j]  # if not convergence

    }
    if(conv1 | conv2) break
  }
  if(i == 1 & !conv1) cw.new = cw
  cw.new = cw.new
  c.new = cw.new * sw
  b.new = sum((zw - Rw %*% cw.new) * sw) / sum(sw)

  return(list(Rw = Rw, z.new = z, zw.new = zw, w.new = w, sw.new = sw, b.new = b.new, c.new = c.new, cw.new = cw.new))
}

sspline.QP = function (R, y, f, lambda0, obj, c.init)
{
  n = length(y)
  mu = obj$linkinv(f)

  # initialize
  w = obj$variance(mu)
  z = f + (y - mu) / w
  b = 0

  zw = z * sqrt(w)
  Rw = R * w
  cw = c.init / sqrt(w)
  sw = sqrt(w)
  cw.new = rep(0, n)
  for(i in 1:10){ # outer iteration
    Dmat = t(R) %*% R + n * lambda0 * R
    dvec = as.vector(t(zw - b * sw) %*% R)
    cw.new = ginv(Dmat) %*% dvec

    loss = abs(cw-cw.new)
    conv = max(loss) < 1e-6

    if(conv) break
    cw = cw.new  # if not convergence
  }

  cw.new = cw.new
  c.new = cw.new * sw
  b.new = sum((zw - Rw %*% cw.new) * sw) / sum(sw)
  return(list(Rw = Rw, z.new = z, zw.new = zw, w.new = w, sw.new = sw, b.new = b.new, c.new = c.new, cw.new = cw.new))
}

# cv.nng = function(model, y, mscale, lambda0, lambda_theta, gamma, obj, algo)
# {
#   cat("-- theta-step -- \n")
#   cat("proceeding... \n")
#   n = length(y)
#   d = length(mscale)
#
#   # solve theta
#   G <- matrix(0, nrow(model$R[, ,1]), d)
#   for (j in 1:d) {
#     G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
#   }
#
#   Gw = G * sqrt(model$w.new)
#   uw = model$zw.new - model$b.new * sqrt(model$w.new) - (n/2) * lambda0 * model$cw.new
#
#   init.theta = rep(1, d)
#
#   len = length(lambda_theta)
#
#   measure <- c(NA, len)
#   save_theta <- list()
#
#
#   for (k in 1:len) {
#     if(algo == "CD") {
#       theta.new = .Call("glm_theta_step", Gw, uw, n, d, init.theta, lambda_theta[k], gamma)
#       # theta.new = nng.cd(Gw, uw, init.theta, lambda_theta[k], gamma)
#       theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
#       save_theta[[k]] <- theta.adj
#     }
#
#     XX = model$zw.new - Gw %*% theta.adj
#     num = t(XX) %*% XX + 1
#     den = (1 - sum(diag( Gw %*% ginv( t(Gw) %*% Gw) %*% t(Gw) )) / n)^2 + 1
#     # measure[k] <- as.vector(num / den /n)
#
#     testfhat = G %*% theta.adj
#     testmu = obj$linkinv(testfhat)
#     measure[k] <- KL(testfhat, obj)
#   }
#
#
#   # plotting error bar
#   if(obj$family == 'gaussian'){
#     main = "Gaussian Family"
#   }
#   if(obj$family == 'binomial'){
#     main = "Binomial Family"
#   }
#   if(obj$family == 'poisson'){
#     main = "Poisson Family"
#   }
#
#   min_id = which.min(measure)
#   optlambda = lambda_theta[min_id]
#
#   ylab = expression("GCV(" * lambda[theta] * ")")
#
#   xrange = log(lambda_theta)
#   plot(xrange, measure, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab, ylim = range(measure), pch = 15, col = 'red')
#
#   theta.new = .Call("glm_theta_step", Gw, uw, n, d, init.theta, optlambda, gamma)
#   # theta.new = nng.cd(Gw, uw, init.theta, optlambda, gamma)
#   # theta.new = save_theta[[id]]
#   theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
#
#   f.new = c(G %*% theta.adj)
#   mu.new = obj$linkinv(f.new)
#
#   if(obj$family == "binomial") miss <- mean(y != ifelse(mu.new < 0.5, 0, 1))
#   if(obj$family == "gaussian") miss <- mean((y - f.new)^2)
#   if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))
#
#   cat("training error:", miss, "\n")
#
#   out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new)
#   return(out)
# }

# model = sspline_cvfit
# lambda0 = model$optlambda
# mscale = wt
cv.nng.subset = function(model, K, y, nbasis, basis.id, mscale, lambda0, lambda_theta, gamma, obj, algo)
{
  cat("-- theta-step -- \n")
  cat("proceeding... \n")
  n = length(y)
  d = length(mscale)

  # solve theta
  G <- matrix(0, n, d)
  for (j in 1:d) {
    G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
  }

  Gw = G * sqrt(model$w.new)
  uw = model$zw.new - model$b.new * sqrt(model$w.new)

  h = rep(0, d)
  for (j in 1:d) {
    h[j] = n * lambda0 * t(model$c.new) %*% model$R[basis.id, , j] %*% model$c.new
  }

  init.theta = rep(1, d)
  len = length(lambda_theta)
  measure <- matrix(NA, 5, len)
  fold = cvsplitID(n, 5, y, family = obj$family)

  # save_theta <- list()
  for(f in 1:5){
    tr_id = as.vector(fold[, -f])
    te_id = fold[, f]

    tr_id = tr_id[!is.na(tr_id)]
    te_id = te_id[!is.na(te_id)]

    tr_n = length(tr_id)
    te_n = length(te_id)

    for (k in 1:len) {
      theta.new = .Call("glm_theta_step", Gw[tr_id,], uw[tr_id], h/2, tr_n, d, init.theta, tr_n * lambda_theta[k] * gamma / 2, tr_n * lambda_theta[k] * (1-gamma))
      theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
      # save_theta[[k]] <- theta.adj

      te_R = array(NA, c(te_n, nbasis, d))
      for(j in 1:d){
        te_R[, , j] = K$K[[j]][te_id, basis.id]
      }

      testfhat = c(wsGram(te_R, theta.adj/mscale^2) %*% model$c.new + model$b.new)

      # testfhat = G[te_id, ] %*% theta.adj
      testmu = obj$linkinv(testfhat)
      # if(obj$family == "gaussian") measure[f, k] <- mean((testfhat - y[te_id])^2)
      # if(obj$family == "binomial") measure[f, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
      # if(obj$family == "poisson") measure[f, k] <- mean(poisson()$dev.resids(y[te_id], testmu, rep(1, te_n)))
      measure[f, k] <- KL(testfhat, obj)
    }
  }
print(measure)
  # plotting error bar
  if(obj$family == 'gaussian'){
    main = "Gaussian Family"
  }
  if(obj$family == 'binomial'){
    main = "Binomial Family"
  }
  if(obj$family == 'poisson'){
    main = "Poisson Family"
  }

  measure_mean = colMeans(measure, na.rm = T)
  measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(te_n)

  sel_id = which(!is.nan(measure_se) & measure_se != Inf)
  measure_mean = measure_mean[sel_id]
  measure_se = measure_se[sel_id]
  lambda_theta = lambda_theta[sel_id]

  min_id = which.min(measure_mean)
  cand_ids = which((measure_mean >= measure_mean[min_id]) &
                     (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
  cand_ids = cand_ids[cand_ids >= min_id]
  std_id = max(cand_ids)
  optlambda = lambda_theta[std_id]

  ylab = expression("GCV(" * lambda[theta] * ")")


  plot(log(lambda_theta), measure_mean, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab,
       ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
  arrows(x0 = log(lambda_theta), y0 = measure_mean - measure_se,
         x1 = log(lambda_theta), y1 = measure_mean + measure_se,
         angle = 90, code = 3, length = 0.1, col = "darkgray")

  theta.new = .Call("glm_theta_step", Gw, uw, h/2, n, d, init.theta, n * optlambda * gamma / 2, n * optlambda * (1-gamma))
  # theta.new = save_theta[[id]]
  theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)

  f.new = c(G %*% theta.adj)
  mu.new = obj$linkinv(f.new)

  if(obj$family == "binomial") miss <- mean(y[te_id] != ifelse(mu.new < 0.5, 0, 1))
  if(obj$family == "gaussian") miss <- mean((y[te_id] - f.new)^2)
  if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, te_n)))

  cat("training error:", miss, "\n")

  out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new)
  return(out)
}

# cv.nng = function(model, K, y, mscale, lambda0, lambda_theta, gamma, obj, algo)
# {
#   cat("-- theta-step -- \n")
#   cat("proceeding... \n")
#   n = length(y)
#   d = length(mscale)
#
#   # solve theta
#   G <- matrix(0, n, d)
#   for (j in 1:d) {
#     G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
#   }
#
#   Gw = G * sqrt(model$w.new)
#   uw = model$zw.new - model$b.new * sqrt(model$w.new)
#
#   h = rep(0, d)
#   for (j in 1:d) {
#     h[j] = n * lambda0 * t(model$c.new) %*% model$R[basis.id, , j] %*% model$c.new
#   }
#
#   init.theta = rep(1, d)
#   len = length(lambda_theta)
#   measure <- matrix(NA, 5, len)
#   fold = cvsplitID(n, 5, y, family = obj$family)
#
#   # save_theta <- list()
#   for(f in 1:5){
#     tr_id = as.vector(fold[, -f])
#     te_id = fold[, f]
#
#     tr_id = tr_id[!is.na(tr_id)]
#     te_id = te_id[!is.na(te_id)]
#
#     tr_n = length(tr_id)
#     te_n = length(te_id)
#
#     for (k in 1:len) {
#       theta.new = .Call("glm_theta_step", Gw[tr_id,], uw[tr_id], h/2, tr_n, d, init.theta, tr_n * lambda_theta[k] * gamma / 2, tr_n * lambda_theta[k] * (1-gamma))
#       theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
#       # save_theta[[k]] <- theta.adj
#
#       te_R = array(NA, c(te_n, nbasis, d))
#       for(j in 1:d){
#         te_R[, , j] = K$K[[j]][te_id, basis.id]
#       }
#
#       testfhat = c(wsGram(te_R, theta.adj/mscale^2) %*% model$c.new + model$b.new)
#
#       # testfhat = G[te_id, ] %*% theta.adj
#       testmu = obj$linkinv(testfhat)
#       # if(obj$family == "gaussian") measure[f, k] <- mean((testfhat - y[te_id])^2)
#       # if(obj$family == "binomial") measure[f, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
#       # if(obj$family == "poisson") measure[f, k] <- mean(poisson()$dev.resids(y[te_id], testmu, rep(1, te_n)))
#       measure[f, k] <- KL(testfhat, obj)
#     }
#   }
#   print(measure)
#   # plotting error bar
#   if(obj$family == 'gaussian'){
#     main = "Gaussian Family"
#   }
#   if(obj$family == 'binomial'){
#     main = "Binomial Family"
#   }
#   if(obj$family == 'poisson'){
#     main = "Poisson Family"
#   }
#
#   measure_mean = colMeans(measure, na.rm = T)
#   measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(te_n)
#
#   sel_id = which(!is.nan(measure_se) & measure_se != Inf)
#   measure_mean = measure_mean[sel_id]
#   measure_se = measure_se[sel_id]
#   lambda_theta = lambda_theta[sel_id]
#
#   min_id = which.min(measure_mean)
#   cand_ids = which((measure_mean >= measure_mean[min_id]) &
#                      (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
#   cand_ids = cand_ids[cand_ids >= min_id]
#   std_id = max(cand_ids)
#   optlambda = lambda_theta[std_id]
#
#   ylab = expression("GCV(" * lambda[theta] * ")")
#
#
#   plot(log(lambda_theta), measure_mean, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab,
#        ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
#   arrows(x0 = log(lambda_theta), y0 = measure_mean - measure_se,
#          x1 = log(lambda_theta), y1 = measure_mean + measure_se,
#          angle = 90, code = 3, length = 0.1, col = "darkgray")
#
#   theta.new = .Call("glm_theta_step", Gw, uw, h/2, n, d, init.theta, n * optlambda * gamma / 2, n * optlambda * (1-gamma))
#   # theta.new = save_theta[[id]]
#   theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
#
#   f.new = c(G %*% theta.adj)
#   mu.new = obj$linkinv(f.new)
#
#   if(obj$family == "binomial") miss <- mean(y[te_id] != ifelse(mu.new < 0.5, 0, 1))
#   if(obj$family == "gaussian") miss <- mean((y[te_id] - f.new)^2)
#   if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, te_n)))
#
#   cat("training error:", miss, "\n")
#
#   out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new)
#   return(out)
# }

nng.cd = function (Gw, uw, theta, lambda_theta, gamma)
{
  n = nrow(Gw)
  d = ncol(Gw)
  r = lambda_theta * gamma * n
  theta.new = rep(0, d)

  for(i in 1:15){
    for(j in 1:d){
      theta.new[j] = 2 * sum((uw - Gw[,-j] %*% theta[-j]) * Gw[,j])

      theta.new[j] = ifelse(theta.new[j] > 0 & r < abs(theta.new[j]), theta.new[j], 0)
      theta.new[j] = theta.new[j] / (sum(Gw[,j]^2) + n * lambda_theta * (1-gamma)) / 2

      loss = abs(theta - theta.new)
      conv = max(loss) < 1e-20

      if(conv) break
      theta[j] = theta.new[j]
    }
    if(conv) break
  }

  if(i == 1 & !conv) theta = rep(0, d)

  return(theta)
}


nng.QP = function (Gw, uw, theta, lambda_theta, gamma)
{
  n = nrow(Gw)
  d = ncol(Gw)
  r = lambda_theta * gamma * n
  theta.new = rep(0, d)

  for(i in 1:15){ # outer iteration
    Dmat = t(Gw) %*% Gw + diag(n * lambda_theta * gamma, d)
    dvec = as.vector(2 * t(uw) %*% Gw)
    Amat = t(rbind(diag(1, d), rep(-1, d)))
    bvec = c(rep(0, d), -lambda_theta)
    theta.new = solve.QP(2 * Dmat, dvec, Amat, bvec)$solution
    theta.new[theta.new < 1e-8] = 0

    loss = abs(theta - theta.new)
    conv = max(loss) < 1e-8

    if(conv) break
    theta = theta.new
  }

  return(theta.new)
}
