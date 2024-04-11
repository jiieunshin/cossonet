############
# mscale = rep(1, d)/wt^2
# cand.lambda = lambda0

cv.sspline = function (x, y, mscale, nfolds, cand.lambda, obj, one.std, type, kparam, algo)
{
  n <- length(y)
  IDmat <- cvsplitID(n, nfolds)

  K = make_anovaKernel(x, x, type = type, kparam)
  d = K$numK
  R = array(NA, c(n, n, d))

  for(j in 1:d){
    R[, , j] = K$K[[j]]
  }

  Rtheta <- wsGram(R, mscale)
  f.init = rep(0.5, n)

  measure <- matrix(NA, ncol = length(cand.lambda), nrow = nfolds)
  miss <- matrix(NA, ncol = length(cand.lambda), nrow = nfolds)
  for (f in 1:nfolds) {
    # print(f)
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    # generate SS-ANOVA
    tr_n = length(trainID)
    te_n = length(testID)

    tr_R = array(NA, c(tr_n, tr_n, d))
    te_R = array(NA, c(te_n, tr_n, d))
    te2_R = array(NA, c(te_n, te_n, d))

    for(j in 1:d){
      tr_R[, , j] = K$K[[j]][trainID, trainID]
      te_R[, , j] = K$K[[j]][testID, trainID]
    }

    tr_Rtheta <- wsGram(tr_R, mscale)
    te_Rtheta <- wsGram(te_R, mscale)
    te2_Rtheta <- wsGram(te2_R, mscale)

    for (k in 1:length(cand.lambda)) {
    # print(k)
      if(algo == "CD"){

        # initialize
        ff = f.init[trainID]
        mu = obj$linkinv(ff)
        w = obj$variance(mu)
        z = ff + (y[trainID] - mu) / w

        c.init = as.vector(glmnet(tr_Rtheta, y[trainID], family = 'gaussian', lambda = cand.lambda[k])$beta)
        zw = z * sqrt(w)
        Rw = tr_Rtheta * w
        cw = c.init / sqrt(w)
        sw = sqrt(w)
        fit = sspline.cd(tr_Rtheta, y[trainID], ff, cand.lambda[k], obj, c.init)

        # fit = .Call("Csspline", zw, Rw, cw, sw, tr_n, cand.lambda[k], PACKAGE = "cdcosso")
        b.new = fit$b.new
        c.new = fit$c.new
        cw.new = fit$cw.new
      }

      if(algo == "QP"){
        c.init = as.vector(glmnet(tr_Rtheta, y[trainID], family = 'gaussian', lambda = cand.lambda[k])$beta)
        fit = sspline.QP(tr_Rtheta, y[trainID], f.init[trainID], cand.lambda[k], obj, c.init)
        b.new = fit$b.new
        c.new = fit$c.new
        cw.new = fit$cw.new
      }
      if(sum(is.nan(cw.new)) == tr_n){
        next
      } else{
        # validation
        testfhat = c(b.new + te_Rtheta %*% c.new)
        testmu = obj$linkinv(testfhat)
        testw = obj$variance(testmu)
        testz = testfhat + (y[testID] - testmu) / testw

        XX = fit$zw.new - fit$Rw %*% fit$cw.new - fit$b.new * fit$sw.new
        num = t(XX) %*% XX
        den = (1 - sum(diag(tr_Rtheta %*% ginv(tr_Rtheta + diag(fit$w.new)/cand.lambda[k]))) / tr_n)^2
        measure[f, k] <- as.vector( num / den / tr_n)

        # measure[f, k] <- mean(KLD(te_y, testfhat, obj))

        if(obj$family == "binomial") miss[f, k] <- mean(ifelse(testmu < 0.5, 0, 1) != y[testID])
        if(obj$family == "gaussian") miss[f, k] <- mean((testmu - y[testID])^2)
        if(obj$family == "poisson") miss[f, k] <- mean(-obj$dev.resids(y[testID], testmu, rep(1, te_n)))
      }
    }
  }

  measure[measure == -Inf | measure == Inf | is.nan(measure)] <- NA
  if(sum(is.na(measure)) == length(cand.lambda) * nfolds){
    conv = FALSE
    return(list(conv = conv))
  }
  rm(tr_Rtheta)
  rm(te_Rtheta)

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

  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22
  cvsd[cvsd == Inf] <- NA
  max_min <- c(min(cvm - cvsd, na.rm = T), max(cvm + cvsd, na.rm = T))
  ylab = expression("GCV(" * lambda[0] * ")")

  # optimal lambda1
  id = which.min(cvm)[1]
  optlambda = cand.lambda[id]

  plot(log(cand.lambda), cvm, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab, ylim = max_min, type = 'n')
  try(arrows(log(cand.lambda), cvm - cvsd, log(cand.lambda), cvm + cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(cand.lambda), cvm, pch = 15, col = 'red')
  abline(v = log(cand.lambda)[id], col = 'darkgrey', lty = 2)

  ###

  miss_cvm <- apply(miss, 2, mean, na.rm = T)
  misS_cvsd <- apply(miss, 2, sd, na.rm = T) / sqrt(nrow(miss)) + 1e-22
  max_min <- c(min(miss_cvm - misS_cvsd, na.rm = TRUE), max(miss_cvm + misS_cvsd, na.rm = TRUE))


  plot(log(cand.lambda), miss_cvm, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = "miss", ylim = max_min, type = 'n')
  try(arrows(log(cand.lambda), miss_cvm - misS_cvsd, log(cand.lambda), miss_cvm + misS_cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(cand.lambda), miss_cvm, pch = 15, col = 'red')
  abline(v = log(cand.lambda)[id], col = 'darkgrey', lty = 2)

  if(algo == "CD"){
    mu = obj$linkinv(f.init)
    w = obj$variance(mu)
    z = f.init + (y - mu) / w

    c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = optlambda)$beta)
    zw = z * sqrt(w)
    Rw = Rtheta * w
    cw = c.init / sqrt(w)
    sw = sqrt(w)

    fit = sspline.cd(Rtheta, y, f.init, optlambda, obj, c.init)

    # fit = .Call("Csspline", zw, Rw, cw, sw, n, optlambda, PACKAGE = "cdcosso")
    f.new = c(fit$b.new + Rtheta %*% fit$c.new)
    mu.new = obj$linkinv(f.new)
    w.new = obj$variance(mu.new)
    z.new = f.new + (y - mu.new) / w.new

    out = list(IDmat = IDmat, measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
               z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = fit$b.new,
               cw.new = fit$cw.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)
  }

  if(algo == "QP"){
    c.init = as.vector(glmnet(Rtheta, y, family = 'gaussian', lambda = optlambda)$beta)
    fit = sspline.QP(Rtheta, y, f.init, optlambda, obj, c.init)
    f.new = c(fit$b.new + Rtheta %*% fit$c.new)
    mu.new = obj$linkinv(f.new)
    w.new = obj$variance(mu.new)
    z.new = f.new + (y - mu.new) / w.new

    out = list(IDmat = IDmat, measure = measure, R = R, w.new = w.new, sw.new = sqrt(w.new),
               z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = fit$b.new,
               cw.new = fit$cw.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)  }

  rm(K)
  rm(Rtheta)

  return(out)
}

# R = tr_R
# y = y[trainID]

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
      conv = max(loss) < 1e-6

      if(conv) break
      cw[j] = cw.new[j]  # if not convergence

    }
    if(conv) break
  }
  if(i == 1 & !conv) cw.new = cw
  cw.new = cw.new
  c.new = cw.new * sw
  b.new = sum((zw - Rw %*% cw.new) * sw) / sum(sw)
  return(list(Rw = Rw, z.new = z, zw.new = zw, w.new = w, sw.new = sw, b.new = b.new, c.new = c.new, cw.new = cw.new))
}

# R = tr_Rtheta
# y = y[trainID]
# f = f.init[trainID]
# lambda0 = cand.lambda[k]
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

# LHS = t(R1) %*% R1 + 2 * n * lambda0 * R2
# RHS = t(R1) %*% z.old
# c.new = ginv(LHS) %*% RHS

# model = sspline_cvfit
# mscale = wt
# lambda0 = sspline_cvfit$optlambda
# init.theta = ifelse(init.theta < 0.5, 0, 1)
# mean(y != ifelse(obj$linkinv(c(G %*% init.theta)) < 0.5, 0, 1))   # 이거는 잘됨. init.theta는 이거로 고정

cv.nng = function(model, x, y, mscale, lambda0, lambda_theta, gamma, nfolds, obj, one.std, algo)
{
  n = length(y)
  d = length(mscale)
  IDmat = model$IDmat

  # solve theta
  G <- matrix(0, nrow(model$R[, ,1]), d)
  for (j in 1:d) {
    G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
  }

  Gw = G * sqrt(model$w.new)
  uw = model$zw.new - model$b.new * sqrt(model$w.new) - (n/2) * lambda0 * model$cw.new
  # init.theta = as.vector(glmnet(Gw, uw, family = "gaussian", lambda = lambda_theta[1])$beta)
  init.theta = rep(1, d)
  len = length(lambda_theta)
  measure <- matrix(0, ncol = len, nrow = nfolds)
  miss <- matrix(0, ncol = len, nrow = nfolds)

  for (f in 1:nfolds) {
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    tr_G = G[trainID,]
    te_G = G[testID,]

    tr_n = length(trainID)
    te_n = length(testID)
    for (k in 1:len) {
      if(algo == "CD") {
        theta.new = nng.cd(Gw[trainID,], uw[trainID], theta = init.theta, lambda_theta[k], gamma)
        # theta.new = .Call("Cnng", Gw[trainID,], uw[trainID], tr_n, d, init.theta, lambda_theta[k], gamma)
      }

      if(algo == "QP") {
        theta.new = nng.QP(Gw[trainID,], uw[trainID], theta = init.theta, lambda_theta[k], gamma)
      }

      testfhat = c(te_G %*% theta.new)
      testmu = obj$linkinv(testfhat)

      XX = model$z.new[trainID] - tr_G %*% theta.new - model$b.new
      num = t(XX) %*% diag(model$w.new[trainID]) %*% XX
      den = (1 - sum(diag( Gw[trainID,] %*% ginv( t(Gw[trainID,]) %*% Gw[trainID,]) %*% t(Gw[trainID,]) )) / tr_n)^2
      measure[f, k] <- as.vector(num / den /tr_n)
      # measure[f, k] <- mean(KLD(te_y, testfhat, obj))

      if(obj$family == "binomial") miss[f, k] <- mean(ifelse(testmu < 0.5, 0, 1) != y[testID])
      if(obj$family == "gaussian") miss[f, k] <- mean((testmu - y[testID])^2)
      if(obj$family == "poisson") miss[f, k] <- mean(-obj$dev.resids(y[testID], testmu, rep(1, te_n)))
    }
  }
  measure[is.nan(measure)] <- NA

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

  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22

  cvm[is.nan(cvm)] <- NA
  cvsd[is.na(cvsd)] <- 0
  id = which.min(cvm)[1]

  miss_cvm <- apply(miss, 2, mean, na.rm = T)
  misS_cvsd <- apply(miss, 2, sd, na.rm = T) / sqrt(nrow(miss)) + 1e-22
  max_min <- c(min(miss_cvm - misS_cvsd, na.rm = TRUE), max(miss_cvm + misS_cvsd, na.rm = TRUE))

  if(one.std){
    st1_err = cvm[id] + cvsd[id] # minimum cv err
    std.id = max(which(cvm[id:len] <= st1_err & cvm[id] <= cvm[id:len]))
    if(is.na(std.id)){
      std.id = id
      optlambda = lambda_theta[std.id]
    } else{
      std.id = ifelse(std.id > id, std.id, id)
      optlambda = lambda_theta[std.id]
    }
  } else{
    optlambda = lambda_theta[id]
  }

  max_min <- c(min(cvm - cvsd, na.rm = TRUE), max(cvm + cvsd, na.rm = TRUE))
  ylab = expression("GCV(" * lambda[theta] * ")")

  xrange = log(lambda_theta)
  plot(xrange, cvm, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab, ylim = max_min, type = 'n')
  arrows(xrange, cvm - cvsd, xrange, cvm + cvsd, angle = 90, code = 3, length = 0.1, col = 'gray')
  points(xrange, cvm, pch = 15, col = 'red')
  abline(v = xrange[id], col = 'darkgrey')
  # text(log(lambda_theta), par("usr")[4], labels = selm, pos = 1)
  if(one.std) abline(v = xrange[std.id], col = 'darkgrey', lty = 2)

  cvm <- apply(miss, 2, mean, na.rm = T)
  cvsd <- apply(miss, 2, sd, na.rm = T) / sqrt(nrow(miss)) + 1e-22
  max_min <- c(min(miss_cvm - misS_cvsd, na.rm = TRUE), max(miss_cvm + misS_cvsd, na.rm = TRUE))

  plot(log(lambda_theta), cvm, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = "miss", ylim = max_min, type = 'n')
  try(arrows(log(lambda_theta), cvm - cvsd, log(lambda_theta), cvm + cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(lambda_theta), cvm, pch = 15, col = 'red')
  abline(v = log(lambda_theta)[id], col = 'darkgrey', lty = 2)

  if(algo == "CD"){
    theta.new = nng.cd(Gw, uw, theta = init.theta, optlambda, gamma)
    # theta.new = .Call("Cnng", Gw, uw, n, d, init.theta, optlambda, gamma)
  }

  if(algo == "QP"){
    theta.new = nng.QP(Gw, uw, theta = init.theta, optlambda, gamma)
  }
  out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new)
  return(out)
}


# zw = model$zw.new[trainID]
# b = model$b.new
# sw = model$sw.new[trainID]
# cw = model$cw.new[trainID]
# w = model$w.new[trainID]
# G = tr_G
# y = y[trainID]

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
    dvec = t(uw) %*% Gw + n * lambda_theta * (1-gamma)
    Amat = diag(1, d)
    bvec = rep(0, d)
    theta.new = solve.QP(2 * Dmat, 2 * dvec, Amat, bvec)$solution
    theta.new[theta.new < 1e-8] = 0

    loss = abs(theta - theta.new)
    conv = max(loss) < 1e-8

    if(conv) break
    theta = theta.new
  }

  return(theta.new)
}


######################
# object = fit3
# testx = te_x
#
# KLD = function(f, y, family = "binomial"){
#   if(family == 'poisson') D = function(f, y) exp(f) - y*f
#   if(family == 'binomial') D = function(f, y) log(exp(1-f)+1) - y*f
#   if(family == 'Cox') D = function(f, y) log(exp(1-f)+1) - y*f
#
#   return(D(f, y))
# }
