############
# mscale = init.theta/wt^2
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

  c.init = as.vector(glmnet(Rtheta, y, family = obj$family, lambda = cand.lambda[1], alpha = 0)$beta)
  c.init = c(scale(c.init))
  f.init = c(Rtheta %*% c.init)

  measure <- matrix(NA, ncol = length(cand.lambda), nrow = nfolds)
  for (f in 1:nfolds) {
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    # generate SS-ANOVA
    tr_n = length(trainID)
    te_n = length(testID)

    tr_R = array(NA, c(tr_n, tr_n, d))
    te_R = array(NA, c(te_n, tr_n, d))

    for(j in 1:d){
      tr_R[, , j] = K$K[[j]][trainID, trainID]
      te_R[, , j] = K$K[[j]][testID, trainID]
    }

    tr_Rtheta <- wsGram(tr_R, mscale)
    te_Rtheta <- wsGram(te_R, mscale)

    for (k in 1:length(cand.lambda)) {

      if(algo == "CD"){

        # initialize
        ff = f.init[trainID]
        mu = obj$linkinv(ff)
        w = obj$variance(mu)
        z = ff + (y[trainID] - mu) / w
        b = 0
        # print(b)

        zw = z * sqrt(w)
        Rw = tr_Rtheta * w
        cw = c.init[trainID] / sqrt(w)
        sw = sqrt(w)
        # sspline_fit = sspline.cd(tr_Rtheta, y[trainID], f.init[trainID], cand.lambda[k], obj, c.init[trainID])

        sspline_fit = .Call("Csspline", zw, Rw, cw, sw, tr_n, cand.lambda[k], PACKAGE = "cdcosso")
        b.new = sspline_fit$b.new
        c.new = sspline_fit$c.new
        cw.new = sspline_fit$cw.new
        # print(b.new)
      }

      if(algo == "QP"){
        sspline_fit = sspline.QP(tr_Rtheta, y[trainID], f.init[trainID], cand.lambda[k], obj, c.init[trainID])
        b.new = sspline_fit$b.new
        c.new = sspline_fit$c.new
        cw.new = sspline_fit$cw.new
      }

      # validation
      # f.new <- fit$f.new
      testfhat = c(b.new + te_Rtheta %*% c.new)
      testmu = obj$linkinv(testfhat)
      testw = obj$variance(testmu)
      testz = testfhat + (y[testID] - testmu) / testw

      testzw = testz * sqrt(testw)
      testRw = te_Rtheta * testw
      rss <- t(testzw - testRw %*% cw.new - b.new * sqrt(testw)) %*% (testzw - testRw %*% cw.new - b.new * sqrt(testw)) + .1
      S = testRw %*% ginv(t(testRw) %*% testRw) %*% t(testRw)
      df = sum(diag(S))
      measure[f, k] <- rss / (1 - df/length(testID) + .1)^2 / length(testID)

      # if(obj$family == "binomial") measure[f, k] <- mean(ifelse(testmu < 0.5, 0, 1) != y[testID])
      # if(obj$family == "gaussian") measure[f, k] <- mean((testmu - y[testID])^2)
      # if(obj$family == "poisson") measure[f, k] <- mean(KLD(testfhat, y[testID]))
      # print(measure)
    }
  }

  rm(tr_Rtheta)
  rm(te_Rtheta)
  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22

  # optimal lambda1
  id = which.min(cvm)[1]
  optlambda = cand.lambda[id]
  # if(one.std){
    # st1_err = cvm[id] + cvsd[id] # minimum cv err
    # std.id = max(which(cvm[1:id] <= st1_err & cvm[1:id] <= cvm[id]))
    # std.id = ifelse(std.id > id, std.id, id)
    # optlambda = cand.lambda[std.id]
  # } else{
  #   optlambda = cand.lambda[id]
  # }
  #
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

  max_min <- c(min(cvm - cvsd, na.rm = T), max(cvm + cvsd, na.rm = T))

  plot(log(cand.lambda), cvm, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = "generalized cross validation", ylim = max_min, type = 'n')
  try(arrows(log(cand.lambda), cvm - cvsd, log(cand.lambda), cvm + cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(cand.lambda), cvm, pch = 15, col = 'red')
  abline(v = log(cand.lambda)[id], col = 'darkgrey', lty = 2)
  # abline(v = log(cand.lambda)[std.id], col = 'darkgrey', lty = 2)

  # c.init = runif(n, -1e-5, 1e-5)

  if(algo == "CD"){
    ff = c(f.init)
    mu = obj$linkinv(ff)
    w = obj$variance(mu)
    z = ff + (y - mu) / w
    # b = mean(z - Rtheta %*% c.init)
    # print(b)

    zw = z * sqrt(w)
    Rw = Rtheta * w
    cw = c.init / sqrt(w)
    sw = sqrt(w)

    # fit = sspline.cd(Rtheta, y, f.init, optlambda, obj, c.init)

    fit = .Call("Csspline", zw, Rw, cw, sw, n, optlambda, PACKAGE = "cdcosso")
    f.new <- c(fit$b.new + Rtheta %*% fit$c.new)
    mu.new = obj$linkinv(f.new)
    w.new = obj$variance(mu.new)
    z.new = f.new + (y - mu.new) / w.new

    out = list(IDmat = IDmat, measure = measure, R = R, f.new = f.new, zw.new = z.new * w.new, b.new = fit$b.new,
               cw.new = fit$cw.new, c.new = fit$c.new, w.new = w.new, optlambda = optlambda)
  }

  if(algo == "QP"){
    fit = sspline.QP(Rtheta, y, f.init, optlambda, obj, c.init)
    f.new <- c(fit$b.new + Rtheta %*% fit$c.new)
    mu.new = obj$linkinv(f.new)
    w.new = obj$variance(mu.new)
    z.new = f.new + (y - mu.new) / w.new

    out = list(IDmat = IDmat, measure = measure, R = R, f.new = f.new, zw.new = z.new * w.new, b.new = fit$b.new,
               cw.new = fit$cw.new, c.new = fit$c.new, w.new = w.new, optlambda = optlambda)
  }

  rm(K)
  rm(Rtheta)

  return(out)
}

# R = tr_R
# y = y[trainID]

sspline.cd = function (R, y, f, lambda0, obj, c.init)
{

  # print(as.vector(glmnet(R, y, family = obj$family, lambda = lambda0, alpha = 0)$beta))
  n = length(y)

  # initialize
  mu = obj$linkinv(f)
  w = obj$variance(mu)
  z = f + (y - mu) / w
  b = 0

  zw = z * sqrt(w)
  Rw = R * w
  cw = c.init / sqrt(w)
  sw = sqrt(w)
  cw.new = rep(0, n)
  for(i in 1:20){ # outer iteration

    for(j in 1:n){
      L = 2 * sum((zw - Rw[,-j] %*% cw[-j] - b * sw) * Rw[,j]) - n * lambda0 * c(Rw[j,-j] %*% cw[-j])
      R = 2 * sum(Rw[,j]^2) + n * lambda0 * Rw[j,j]
      cw.new[j] = L/R

      loss = abs(cw-cw.new)
      conv = max(loss) < 1e-6

      if(conv) break
      cw = cw.new  # if not convergence

    }
  }
  if(i == 1 & !conv) cw.new = cw
  cw.new = cw.new / sd(cw.new)
  c.new = cw.new * sqrt(w)
  b.new = sum((zw - Rw %*% cw.new) * sw) / sum(sw)

  return(list(w.new = w, b.new = b.new, c.new = c.new, zw.new = z * w, sw.new = sqrt(w), cw.new = cw.new))
}

# sspline.cd = function (R, y, f, lambda0, obj, c.init)
# {
#   n = dim(R)[1]
#
#   # initialize
#   mu = obj$linkinv(f)
#   w = obj$variance(mu)
#   z = f + (y - mu) / w
#   b = 0
#
#   zw = z * sqrt(w)
#   Rw = R * w
#   cw = c.init / sqrt(w)
#   sw = sqrt(w)
#
#   for(i in 1:20){
#     cw.new = sapply(1:n, function(j){
#       L = 2 * mean((zw - Rw[,-j] %*% cw[-j] - b * sw) * Rw[,j]) - n * lambda0 * c(Rw[j,-j] %*% cw[-j])
#       R = 2 * sum(Rw[,j]^2) + n * lambda0 * Rw[j,j]
#       L/R
#     })
#
#     loss = abs(cw-cw.new)
#     conv = max(loss) < 1e-5
#
#     if(conv) break
#     cw = cw.new  # if not convergence
#   }
#   # cw.new = cw.new
#   cw.new = cw.new
#   c.new = cw.new * sqrt(w)
#   b.new = sum((zw - Rw %*% cw.new) * sw) / sum(sw)
#
#   return(list(w.new = w, b.new = b.new, c.new = c.new, zw.new = z * w, sw.new = sqrt(w), cw.new = cw.new))
# }

sspline.QP = function (R, y, f, lambda0, obj, c.init)
{
  n = dim(R)[1]

  # initialize
  mu = obj$linkinv(f)
  w = obj$variance(mu)
  z = f + (y - mu) / w
  b = mean(z - R %*% c.init)

  zw = z * sqrt(w)
  Rw = R * w
  cw = c.init / sqrt(w)
  sw = sqrt(w)

  # iteration
  cw = c.init
  for(i in 1:20){

    D = (t(Rw) %*% Rw + n * lambda0 * Rw)
    dvec = t(Rw) %*% (zw - b*sw)
    cw.new = MASS::ginv(D) %*% dvec

    loss = abs(cw-cw.new)
    conv = max(loss) < 1e-5

    if(conv) break
    cw =  cw.new  # if not convergence
  }
  # print(conv)
  # cw.new = c(scale(cw.new))
  c.new = c(cw.new * sqrt(w))
  b.new = sum((zw - Rw %*% cw.new) * sw) / sum(sw)
  # cat("convergence", conv, "iteration", i, "\n")
  return(list(w.new = w, b.new = b.new, c.new = c.new, zw.new = z * w, sw.new = sqrt(w), cw.new = c(cw.new)))
}

# LHS = t(R1) %*% R1 + 2 * n * lambda0 * R2
# RHS = t(R1) %*% z.old
# c.new = ginv(LHS) %*% RHS

# model = sspline_cvfit
# mscale = wt
# lambda0 = optlambda0
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
  init.theta = as.vector(glmnet(Gw, uw, family = "gaussian", lambda = lambda0)$beta)
  len = length(lambda_theta)
  measure <- matrix(0, ncol = len, nrow = nfolds)
  l = 0
  for (f in 1:nfolds) {
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    tr_G = G[trainID,]
    te_G = G[testID,]

    tr_n = length(trainID)
    te_n = length(testID)
    for (k in 1:len) {
      l = l + 1
      if(algo == "CD") {
        # theta.new = nng.cd(Gw[trainID,], uw[trainID], theta = init.theta, lambda_theta[k], gamma)

        theta.new = .Call("Cnng", Gw[trainID,], uw[trainID], tr_n, d, init.theta, lambda_theta[k], gamma)
        # print(theta.new)
      }

      if(algo == "QP") {
        theta.new = nng.QP(model$zw.new[trainID], model$b.new, model$cw.new[trainID], model$w.new[trainID], tr_G,
                         theta = init.theta, lambda0, lambda_theta[k], gamma)
      }

      testfhat = c(te_G %*% theta.new)
      testmu = obj$linkinv(testfhat)
      testw = obj$variance(testmu)
      testz = testfhat + (y[testID] - testmu) / testw
      testzw = testz * sqrt(testw)
      testGw = te_G * sqrt(testw)
      testuw = testzw - model$b.new * sqrt(testw) - (te_n/2) * lambda0 * model$cw.new[testID]
      rss <- t(testuw - testGw %*% theta.new) %*% (testuw - testGw %*% theta.new) + .1
      l1 = gamma * sum(abs(theta.new)) + (1-gamma) * norm(theta.new, "2")
      l2 = gamma * sum(abs(ginv(theta.new))) + (1-gamma) * norm(ginv(theta.new), "2")
      S = l1 + l2
      measure[f, k] <- rss / (1 - d * S/te_n + .1)^2 / te_n

      # if(obj$family == "binomial") measure[f, k] <- mean(ifelse(testmu < 0.5, 0, 1) != y[testID])
      # if(obj$family == "gaussian") measure[f, k] <- mean((testmu - y[testID])^2)
      # if(obj$family == "poisson") measure[f, k] <- mean(KLD(testfhat, y[testID]))
    }
  }
  measure[is.nan(measure)] <- NA

  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22

  cvm[is.nan(cvm)] <- NA
  cvsd[is.na(cvsd)] <- 0
  # selm = floor(apply(sel, 2, mean))
  id = which.min(cvm)[1]

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
  max_min <- c(min(cvm - cvsd, na.rm = TRUE), max(cvm + cvsd, na.rm = TRUE))

  xrange = log(lambda_theta)
  plot(xrange, cvm, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = "generalized cross validation", ylim = max_min, type = 'n')
  arrows(xrange, cvm - cvsd, xrange, cvm + cvsd, angle = 90, code = 3, length = 0.1, col = 'gray')
  points(xrange, cvm, pch = 15, col = 'red')
  abline(v = xrange[id], col = 'darkgrey')
  # text(log(lambda_theta), par("usr")[4], labels = selm, pos = 1)
  if(one.std) abline(v = xrange[std.id], col = 'darkgrey')

  if(algo == "CD"){
    # theta.new = nng.cd(Gw, uw, theta = init.theta, optlambda, gamma)
    theta.new = .Call("Cnng", Gw, uw, n, d, init.theta, optlambda, gamma)
    f.new = c(G %*% as.matrix(theta.new))
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new, f.new = f.new)
  }

  if(algo == "QP"){
    theta.new = nng.QP(model$zw.new, model$b.new, model$cw.new, model$w.new, G,
                     init.theta, lambda0, optlambda, gamma, obj)
    f.new = c(G %*% as.matrix(theta.new))
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new, f.new = f.new)
  }

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

  for(i in 1:20){

    for(j in 1:d){
      theta.new[j] = 2 * sum((uw - Gw[,-j] %*% theta[-j]) * Gw[,j])

      theta.new[j] = ifelse(theta.new[j] > 0 & r < abs(theta.new[j]), theta.new[j], 0)
      theta.new[j] = theta.new[j] / (sum(Gw[,j]^2) + n * lambda_theta * (1-gamma)) / 2

      loss = abs(theta[j] - theta.new[j])
      conv = max(loss) < 1e-20

      if(i != 1 & conv) break
      theta = theta.new
    }
  }

  if(i == 1 & !conv){
    theta = rep(0, d)
  }
  # else if(sum(theta.new = 0) == d){
  #   theta = theta.new / sd(theta.new)
  # }
  # print(theta)
  # if(sum(theta == 0) != d) theta = theta / sd(theta)

  return(theta)
}


nng.QP = function (zw, b, cw, w, G, theta, lambda0, lambda_theta, gamma, obj)
{
  n = nrow(G)
  d = ncol(G)

  Gw = G * sqrt(w)
  uw = zw - b * sqrt(w) - (n/2) * lambda0 * cw

  # print(algo)
  Amat = diag(1, d)
  bvec = rep(0, d)
  for(i in 1:20){
    Dmat = 2 * (t(Gw) %*% Gw  + diag(lambda_theta * (1-gamma), d))
    dvec = c(2 * t(uw) %*% Gw - gamma * lambda_theta)

    # dvec = ifelse(dvec > 0 & abs(dvec) > lambda_theta * gamma, dvec, 0)
    # theta.new = c(ginv(Dmat) %*% dvec)
    theta.new <- solve.QP(Dmat, dvec, t(Amat), bvec)$solution

    theta.new[theta.new < 1e-10] <- 0
    if(sum(theta.new == 0) < d) theta.new = theta.new / sd(theta.new)

    loss = abs(theta-theta.new)
    conv = max(loss) < 1e-10
    if(conv) break

    theta = theta.new
  }
  return(theta.new)
}


######################
# object = fit3
# testx = te_x

KLD = function(f, y, family = "binomial"){
  if(family == 'poisson') D = function(f, y) exp(f) - y*f
  if(family == 'binomial') D = function(f, y) log(exp(1-f)+1) - y*f
  if(family == 'Cox') D = function(f, y) log(exp(1-f)+1) - y*f

  return(D(f, y))
}
