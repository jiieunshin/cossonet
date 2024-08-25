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
#       fit = .Call("glm_c_step", zw, Rw, cw, sw, n, cand.lambda[k], PACKAGE = "cdcosso")
#       b.new = fit$b.new
#       c.new = fit$c.new
#       cw.new = fit$cw.new
#     }
#
#     if(algo == "QP"){
#       c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = cand.lambda[k])$beta)
#       fit = sspline.QP(Rtheta, y, f.init, cand.lambda[k], obj, c.init)
#       b.new = fit$b.new
#       c.new = fit$c.new
#       cw.new = fit$cw.new
#     }
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
#       XX = fit$zw.new - Rw %*% fit$cw.new - fit$b.new * sqrt(w)
#       num = t(XX) %*% XX + 1
#       # den = (1 - sum(diag(Rtheta %*% ginv(Rtheta + diag(w)/cand.lambda[k]))) / n)^2
#       S = Rw %*% ginv(t(Rw) %*% Rw) %*% t(Rw)
#       den = (1 - sum(diag(S)) / n)^2 + 1
#       measure[k] <- as.vector( num / den / n)
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
#   if(algo == "CD"){
#     mu = obj$linkinv(f.init)
#     w = obj$variance(mu)
#     z = f.init + (y - mu) / w
#
#     c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = optlambda)$beta)
#     zw = z * sqrt(w)
#     Rw = Rtheta * w
#     cw = c.init / sqrt(w)
#     sw = sqrt(w)
#
#     fit = .Call("glm_c_step", zw, Rw, cw, sw, n, optlambda, PACKAGE = "cdcosso")
#     f.new = c(fit$b.new + Rtheta %*% fit$c.new)
#     mu.new = obj$linkinv(f.new)
#     w.new = obj$variance(mu.new)
#     z.new = f.new + (y - mu.new) / w.new
#
#     if(obj$family == "binomial") miss <- mean(y != ifelse(mu.new < 0.5, 0, 1))
#     if(obj$family == "gaussian") miss <- mean((y - f.new)^2)
#     if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))
#
#     cat("training error:", miss, "\n")
#
#     out = list(measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
#                z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = fit$b.new,
#                cw.new = fit$cw.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)
#   }
#
#   if(algo == "QP"){
#     c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = optlambda)$beta)
#     fit = sspline.QP(Rtheta, y, f.init, optlambda, obj, c.init)
#     f.new = c(fit$b.new + Rtheta %*% fit$c.new)
#     mu.new = obj$linkinv(f.new)
#     w.new = obj$variance(mu.new)
#     z.new = f.new + (y - mu.new) / w.new
#
#     out = list(measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
#                z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = fit$b.new,
#                cw.new = fit$cw.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)  }
#
#   rm(K)
#   rm(Rtheta)
#
#   return(out)
# }

# cand.lambda = lambda0
# mscale = wt
# obj = gaussian()
cv.sspline = function (K, y, mscale, cand.lambda, obj, type, kparam, algo, show)
{
  cat("-- c-step -- \n")
  cat("proceeding... \n")
  d = K$numK
  n <- length(y)
  len = length(cand.lambda)

  R = array(NA, c(n, n, d))
  for(j in 1:d){
    R[, , j] = K$K[[j]]
  }

  Rtheta <- combine_kernel(R, mscale)

  # initialize
  f.init = rep(0.5, n)
  ff = f.init
  mu = obj$linkinv(ff)
  w = obj$variance(mu)
  z = ff + (y - mu) / w

  #
  fold = cvsplitID(n, 5)
  measure <- matrix(NA, 5, len)
  for(f in 1:5){
    if(obj$family != "binomial"){
      tr_id = as.vector(fold[, -f])
      te_id = fold[, f]
      m = length(tr_id)
    } else{

    }

    tr_R = array(NA, c(m, n, d))
    for(j in 1:d){
      tr_R[, , j] = K$K[[j]][tr_id, ]
    }

    tr_Rtheta <- combine_kernel(tr_R, mscale)

    te_R = array(NA, c(n-m, n, d))
    for(j in 1:d){
      te_R[, , j] = K$K[[j]][te_id, ]
    }

    te_Rtheta <- combine_kernel(te_R, mscale)

    #
    zw = z[tr_id] * sqrt(w[tr_id])
    Rw = tr_Rtheta * w[tr_id]
    sw = sqrt(w)[tr_id]
    Rw2 = Rtheta * w
    for (k in 1:len){

      if(algo == "CD"){
        c.init = as.vector(glmnet(tr_Rtheta, y[tr_id], family = 'gaussian', lambda = 1e-4)$beta)

        cw = c.init / sqrt(w)

        fit = .Call("glm_c_step", zw, Rw, Rw2, cw, sw, m, n, cand.lambda[k], PACKAGE = "cdcosso")
        b.new = fit$b.new
        cw.new = fit$cw.new
        c.new = cw.new * sqrt(w)
        # cat("R calculate:", sum(zw - Rw %*% cw.new) / sum(sw), "\n")
        # cat("C calculate:", b.new, "\n")
      }

      # if(algo == "QP"){
      #   c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = cand.lambda[k])$beta)
      #   fit = sspline.QP(Rtheta, y, f.init, cand.lambda[k], obj, c.init)
      #   b.new = fit$b.new
      #   c.new = fit$c.new
      #   cw.new = fit$cw.new
      # }
      if(sum(is.nan(cw.new)) == n){
        next
      } else{

        # validation
        testfhat = c(b.new + te_Rtheta %*% c.new)
        testmu = obj$linkinv(testfhat)

        # XX = zw - Rw %*% cw.new - b.new * sqrt(w[tr_id])
        # num = t(XX) %*% XX + 1
        # den = (1 - sum(diag(Rw %*% ginv(Rw + diag(w[tr_id])/cand.lambda[k]))) / n)^2
        # S = Rw %*% ginv(t(Rw) %*% Rw) %*% t(Rw)
        # den = (1 - sum(diag(S)) / n)^2 + 1
        # measure[f, k] <- num / den
        if(obj$family == "gaussian") measure[f, k] <- mean((testfhat - y[te_id])^2)
        # if(obj$family == "binomial") measure[f, k] <- mean(y[te_id] != ifelse(mu.new < 0.5, 0, 1))
        # if(obj$family = "poisson") measure[k] <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))

      }
    }
  }
# print(measure)
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
  measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(n-m)
  min_id = which.min(measure_mean)
  cand_ids = which((measure_mean >= measure_mean[min_id]) &
                     (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
  cand_ids = cand_ids[cand_ids >= min_id]
  std_id = max(cand_ids)
  optlambda = cand.lambda[std_id]

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

  if(algo == "CD"){
    # f.init = rep(0.5, n)
    # ff = f.init
    # mu = obj$linkinv(ff)
    # w = obj$variance(mu)
    # z = ff + (y - mu) / w

    c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = optlambda)$beta)

    cw = c.init / sqrt(w)

    zw = z * sqrt(w)
    Rw = Rtheta * w
    sw = sqrt(w)

    fit = .Call("glm_c_step", zw, Rw, Rw, cw, sw, n, n, optlambda, PACKAGE = "cdcosso")
    b.new = fit$b.new
    cw.new = fit$cw.new
    c.new = fit$cw.new * sqrt(w)

    f.new = c(b.new + Rtheta %*% c.new)
    mu.new = obj$linkinv(f.new)
    w.new = obj$variance(mu.new)
    z.new = f.new + (y - mu.new) / w.new

    if(obj$family == "binomial") miss <- mean(y != ifelse(mu.new < 0.5, 0, 1))
    if(obj$family == "gaussian") miss <- mean((y - f.new)^2)
    if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))

    cat("training error:", miss, "\n")

    out = list(measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
               z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = b.new,
               cw.new = cw.new, c.new = c.new, optlambda = optlambda, conv = TRUE)
  }

  if(algo == "QP"){
    c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = optlambda)$beta)
    fit = sspline.QP(Rtheta, y, f.init, optlambda, obj, c.init)
    f.new = c(fit$b.new + Rtheta %*% fit$c.new)
    mu.new = obj$linkinv(f.new)
    w.new = obj$variance(mu.new)
    z.new = f.new + (y - mu.new) / w.new

    out = list(measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
               z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = b.new,
               cw.new = cw.new, c.new = c.new, optlambda = optlambda, conv = TRUE)
    }

  rm(K)
  rm(Rtheta)

  return(out)
}


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
  for(i in 1:10){ # outer iteration

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

cv.nng = function(model, y, mscale, lambda0, lambda_theta, gamma, obj, algo)
{
  cat("-- theta-step -- \n")
  cat("proceeding... \n")
  n = length(y)
  d = length(mscale)

  # solve theta
  G <- matrix(0, nrow(model$R[, ,1]), d)
  for (j in 1:d) {
    G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
  }

  Gw = G * sqrt(model$w.new)
  uw = model$zw.new - model$b.new * sqrt(model$w.new) - (n/2) * lambda0 * model$cw.new

  init.theta = rep(1, d)

  if(algo == "QP") lambda_theta = exp(seq(log(1e-4), log(40), length.out = length(lambda_theta)))
  len = length(lambda_theta)

  measure <- matrix(NA, 5, len)
  save_theta <- list()
  for(f in 1:5){
    if(obj$family != "binomial"){
      fold = cvsplitID(n, 5)
      tr_id = as.vector(fold[, -f])
      te_id = fold[, f]
      m = length(tr_id)
    } else{
    }

    for (k in 1:len) {
      if(algo == "CD") {
        theta.new = .Call("glm_theta_step", Gw[tr_id,], uw[tr_id], m, d, init.theta, lambda_theta[k], gamma)
        theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
        save_theta[[k]] <- theta.adj
      }

      # if(algo == "QP") {
      #   theta.new = nng.QP(Gw, uw, theta = init.theta, lambda_theta[k], gamma)
      #   save_theta[[k]] <- theta.new
      # }

      # XX = model$zw.new - Gw %*% theta.adj
      # num = t(XX) %*% XX + 1
      # den = (1 - sum(diag( Gw %*% ginv( t(Gw) %*% Gw) %*% t(Gw) )) / n)^2 + 1
      # measure[k] <- as.vector(num / den /n)

      testfhat = G[te_id, ] %*% theta.adj
      testmu = obj$linkinv(testfhat)
      if(obj$family == "gaussian") measure[f, k] <- mean((testfhat - y[te_id])^2)
      if(obj$family == "binomial") measure[f, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
      # if(obj$family = "poisson") measure[k] <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))
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
# print(save_theta)
  # min_val = min(measure)
  # id = max(which(measure == min_val))
  # optlambda = lambda_theta[id]

  measure_mean = colMeans(measure, na.rm = T)
  measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(n-m)
  min_id = which.min(measure_mean)
  cand_ids = which((measure_mean >= measure_mean[min_id]) &
                     (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
  cand_ids = cand_ids[cand_ids >= min_id]
  std_id = max(cand_ids)
  optlambda = lambda_theta[std_id]

  ylab = expression("GCV(" * lambda[theta] * ")")

  xrange = log(lambda_theta)
  # plot(xrange, measure, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab, ylim = range(measure), pch = 15, col = 'red')

  plot(log(lambda_theta), measure_mean, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab,
       ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
  arrows(x0 = log(lambda_theta), y0 = measure_mean - measure_se,
         x1 = log(lambda_theta), y1 = measure_mean + measure_se,
         angle = 90, code = 3, length = 0.1, col = "darkgray")

  if(algo == "CD"){
    theta.new = .Call("glm_theta_step", Gw, uw, n, d, init.theta, optlambda, gamma)
    # theta.new = save_theta[[id]]
    theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)

    f.new = c(G %*% theta.adj)
    mu.new = obj$linkinv(f.new)

    if(obj$family == "binomial") miss <- mean(y != ifelse(mu.new < 0.5, 0, 1))
    if(obj$family == "gaussian") miss <- mean((y - f.new)^2)
    if(obj$family == "poisson") miss <- mean(poisson()$dev.resids(y, mu.new, rep(1, n)))

    cat("training error:", miss, "\n")
  }

  if(algo == "QP"){
    theta.new = nng.QP(Gw, uw, theta = init.theta, optlambda, gamma)
    # theta.new = save_theta[[id]]
  }
  out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new)
  return(out)
}

nng.cd = function (Gw, uw, theta, lambda_theta, gamma)
{
  n = nrow(Gw)
  d = ncol(Gw)
  r = lambda_theta * gamma * n
  theta.new = rep(0, d)

  for(i in 1:10){
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

  for(i in 1:10){ # outer iteration
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
