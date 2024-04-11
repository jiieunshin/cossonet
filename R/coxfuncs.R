RiskSet = function (time, status)
{
  uniqTime = sort(unique(time[status == 1]))
  RiskSet = matrix(0, ncol = length(uniqTime), nrow = length(time))
  for (k in 1:length(uniqTime)) {
    risk.id = which(time >= uniqTime[k])
    RiskSet[risk.id, k] = risk.id
  }
  return(RiskSet)
}

PartialLik = function (time, status, RS, K, a, neg = FALSE) {
  pl = rep(NA, ncol(RS))
  eventtime = unique(time[status == 1])
  tie.size = as.numeric(table(time[status == 1]))
  for (k in 1:ncol(RS)) {
    failid = which(time == eventtime[k])
    pl[k] = tie.size[k] * log(sum(exp(K[RS[,  k],] %*% a), na.rm = T) + 1e-10)
  }
  pl = sum(pl) - t(status) %*% K %*% a

  if(neg) pl = -pl
  return(pl)
}

# time = unlist(y[, 'time'])
# status = unlist(y[, 'status'])
# mscale = rep(1, d)/wt^2
# nfolds = 5
# cand.lambda = lambda0
cv.getc = function(x, time, status, mscale, nfolds, cand.lambda, one.std, type, kparam, algo)
{
  n <- length(time)
  IDmat <- cvsplitID(n, nfolds)

  K = make_anovaKernel(x, x, type = type, kparam)
  d = K$numK
  R = array(NA, c(n, n, d))

  for(j in 1:d){
    R[, , j] = K$K[[j]]
  }

  Rtheta <- wsGram(R, mscale)

  RS = RiskSet(time, status)

  measure <- matrix(NA, ncol = length(cand.lambda), nrow = nfolds)
  miss <- matrix(NA, ncol = length(cand.lambda), nrow = nfolds)
  for (f in 1:nfolds) {
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    # generate SS-ANOVA
    tr_n = length(trainID)
    te_n = length(testID)

    tr_RS = RiskSet(time[trainID], status[trainID])
    te_RS = RiskSet(time[testID], status[testID])

    tr_R = array(NA, c(tr_n, tr_n, d))
    te_R = array(NA, c(te_n, tr_n, d))

    for(j in 1:d){
      tr_R[, , j] = K$K[[j]][trainID, trainID]
      te_R[, , j] = K$K[[j]][testID, trainID]
    }

    tr_Rtheta <- wsGram(tr_R, mscale)
    te_Rtheta <- wsGram(te_R, mscale)

    for (k in 1:length(cand.lambda)){
      # dyn.load("src/coxfuncs.dll")
      # .Call("Cget_c", tr_Rtheta, Rtheta, n, tr_n, tr_RS, c.init, cand.lambda[k])
      if(algo == "CD"){
        c.init = as.vector(glmnet(tr_Rtheta, cbind(time = time[trainID], status = status[trainID]), family = 'cox', lambda = cand.lambda[k], alpha = 0)$beta)
        fit = getc.cd(tr_Rtheta, c.init, time[trainID], status[trainID], cand.lambda[k], tr_RS)
      }

      if(algo == "QP"){
        # fit = getc.QP(tr_Rtheta, Rtheta, time[trainID], status[trainID], tr_RS, cand.lambda[k])
      }

      Lik = PartialLik(time[trainID], status[trainID], tr_RS, tr_Rtheta, fit$c.new, neg = FALSE)

      XX = fit$zw.new - fit$Rw %*% fit$cw.new - fit$b.new * fit$sw.new
      num = t(XX) %*% XX
      den = (1 - sum(diag(tr_Rtheta %*% ginv(tr_Rtheta + diag(fit$w.new)/cand.lambda[k]))) / tr_n)^2
      measure[f, k] <- as.vector(num / den / tr_n)

      miss[f, k] = Lik
    }
  }
  # print(measure)
  rm(tr_Rtheta)
  rm(te_Rtheta)
  measure[measure == -Inf | measure == Inf | is.nan(measure)] <- NA
  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22

  # optimal lambda1
  id = which.min(cvm)[1]
  optlambda = cand.lambda[id]

  # plotting error bar
  main = "Cox family"
  max_min <- c(min(cvm - cvsd), max(cvm + cvsd))
  ylab = expression("GCV(" * lambda[0] * ")")

  plot(log(cand.lambda), cvm, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab, ylim = max_min, type = 'n')
  try(arrows(log(cand.lambda), cvm - cvsd, log(cand.lambda), cvm + cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(cand.lambda), cvm, pch = 15, col = 'red')
  abline(v = log(cand.lambda)[id], col = 'darkgrey', lty = 2)
  # if(one.std) abline(v = log(cand.lambda)[std.id], col = 'darkgrey')

  ###
  miss_cvm <- apply(miss, 2, mean, na.rm = T)
  misS_cvsd <- apply(miss, 2, sd, na.rm = T) / sqrt(nrow(miss)) + 1e-22
  max_min <- c(min(miss_cvm - misS_cvsd, na.rm = TRUE), max(miss_cvm + misS_cvsd, na.rm = TRUE))

  plot(log(cand.lambda), miss_cvm, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = "partial likelihood", ylim = max_min, type = 'n')
  try(arrows(log(cand.lambda), miss_cvm - misS_cvsd, log(cand.lambda), miss_cvm + misS_cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(cand.lambda), miss_cvm, pch = 15, col = 'red')
  abline(v = log(cand.lambda)[id], col = 'darkgrey', lty = 2)
  ###

  c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox', lambda = optlambda, alpha = 0)$beta)
  fit = getc.cd(Rtheta, c.init, time, status, optlambda, RS)
  # f.new = c(Rtheta %*% c.new)
  out = list(IDmat = IDmat, RS = RS, measure = measure, R = R, w.new = fit$w.new, zw.new = fit$zw.new, b.new = fit$b.new,
             cw.new = fit$cw.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)

  rm(K)
  rm(Rtheta)

  return(out)
}


# Rtheta = tr_Rtheta
# time = time[trainID]
# status = status[trainID]
# lambda0 = cand.lambda[k]
# Risk = tr_RS
getc.cd = function(Rtheta, c.init, time, status, lambda0, Risk)
{
  n = ncol(Rtheta)
  wz = calculate_wz_for_c(c.init, Rtheta, time, status, Risk)
  w = wz$weight
  z = wz$z
  b = 0
  rm(wz)

  zw = z * sqrt(w)
  Rw = Rtheta * w
  cw = c.init
  cw.new = temp = c.init / sqrt(w)
  sw = sqrt(w)

  for(i in 1:20){ # outer iteration
    for(j in 1:n){
      L = 2 * sum((zw - Rw[,-j] %*% cw[-j] - b * sw) * Rw[,j]) - n * lambda0 * c(Rw[j,-j] %*% cw[-j])
      R = 2 * sum(Rw[,j]^2) + n * lambda0 * Rw[j,j]
      temp[j] = L/R

      loss = abs(cw - temp)
      conv1 = max(loss) < 1e-6

      if(conv1) break
      cw[j] <- cw.new[j] <- temp[j]
    }
    if(conv1) break
  }

  c.new = cw.new * sw
  b.new = sum((zw - Rw %*% cw.new) * sw) / sum(sw)

  return(list(Rw = Rw, zw.new = zw, w.new = w, sw.new = sw, b.new = b.new, c.new = c.new, cw.new = cw.new))
}

# R = tr_Rtheta
# time = time[trainID]
# status = status[trainID]
# lambda0 = lambda0[1]
# Risk = tr_RS
# calculate_GH_for_c = function (c.init, R, time, status, lambda0, Risk)
# {
#   n = length(time)
#   tie.size = as.numeric(table(time[status == 1]))
#   # if (min(eigen(R)$value) < 0)
#   #   R = R + 1e-08 * diag(nrow(R))
#   eta = R %*% c.init
#   Hess.FullNumer.unScale = array(NA, dim = c(length(c.init), length(c.init), n))
#   for (i in 1:n) Hess.FullNumer.unScale[, , i] = R[i, ] %*% t(R[i, ])
#
#   Grad.Term1 = -t(R) %*% status/n
#   Grad.Term2 = matrix(NA, ncol = ncol(Risk), nrow = length(c.init))
#   Grad.Term3 = 2 * lambda0 * R %*% c.init
#   Grad.FullNumer = t(R) %*% diag(as.numeric(exp(eta)))
#   Grad.FullDenom = Hess.FullDenom = exp(eta)
#   Hess.FullNumer = Hess.FullNumer.unScale * array(rep(exp(eta), each = length(c.init)^2),
#                                                   dim = c(length(c.init), length(c.init), n)
#                                                   )
#   Hess.Term1 = Hess.Term2 = array(NA, dim = c(length(c.init), length(c.init), ncol(Risk)))
#   k = 1
#   tempSum.exp.eta = sum(exp(eta[Risk[, k]]), na.rm = TRUE)
#   temp.Gradient.numer = apply(Grad.FullNumer[, Risk[, k]], 1, sum, na.rm = TRUE)
#   temp.Hessian.numer = apply(Hess.FullNumer[, , Risk[, k]], c(1, 2), sum, na.rm = TRUE)
#   Grad.Term2[, k] = tie.size[k] * temp.Gradient.numer/tempSum.exp.eta
#   Hess.Term1[, , k] = temp.Hessian.numer/tempSum.exp.eta
#   Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
#   for (k in 2:ncol(Risk)) {
#     excludeID = Risk[, k - 1][!Risk[, k - 1] %in% Risk[, k]]
#     tempSum.exp.eta = tempSum.exp.eta - sum(exp(eta[excludeID]))
#     if (length(excludeID) > 1) {
#       temp.Gradient.numer = temp.Gradient.numer - apply(Grad.FullNumer[, excludeID], 1, sum)
#       temp.Hessian.numer = temp.Hessian.numer - apply(Hess.FullNumer[, , excludeID], c(1, 2), sum)
#     }
#     else {
#       temp.Gradient.numer = temp.Gradient.numer - Grad.FullNumer[, excludeID]
#       temp.Hessian.numer = temp.Hessian.numer - Hess.FullNumer[, , excludeID]
#     }
#     Grad.Term2[, k] = tie.size[k] * temp.Gradient.numer/tempSum.exp.eta
#     Hess.Term1[, , k] = temp.Hessian.numer/tempSum.exp.eta
#     Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
#   }
#   Grad.Term2 = apply(Grad.Term2, 1, sum)/n
#   Gradient = as.vector(Grad.Term1 + Grad.Term2 + Grad.Term3)
#   Hessian = apply(Hess.Term1, c(1, 2), sum)/n - apply(Hess.Term2, c(1, 2), sum)/n + 2 * lambda0 * R
#   return(list(Gradient = Gradient, Hessian = Hessian))
# }

calculate_wz_for_c = function(c.init, R, time, status, RS){
  n = length(time)
  weight = z = rep(0, n)

  for (k in 1:n) {
    Sum.exp.eta.Grad = Sum.exp.eta.Hess = 0
    id = which(RS[k,] > 0)
    eta = as.numeric(R[k,] %*% c.init)
    exp.eta = exp(eta)

    for(r in id){
      Sum.exp.eta = sum(exp(R[RS[,r],] %*% c.init))
      Sum.exp.eta.Grad = Sum.exp.eta.Grad + exp.eta / Sum.exp.eta # {j in R_i} exp(R_j c)
      Sum.exp.eta.Hess = Sum.exp.eta.Hess + ( exp.eta * Sum.exp.eta - exp.eta^2 ) / Sum.exp.eta^2
    }

    Grad.Term = status[k] - Sum.exp.eta.Grad
    weight[k] = Sum.exp.eta.Hess
    z[k] = eta + Grad.Term / weight[k]
  }

  return(list(z = z, weight = weight))
}

# model = getc_cvfit
# lambda0 = getc_cvfit$optlambda
# mscale = wt
cv.gettheta = function (model, x, time, status, mscale, lambda0, lambda_theta, gamma, nfolds, one.std, type, kparam, algo){
  n = length(time)
  d = length(mscale)
  IDmat = model$IDmat

  # solve theta
  G <- matrix(0, nrow(model$R[, ,1]), d)
  for (j in 1:d) {
    G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
  }

  # Gw = G * sqrt(model$w.new)
  # uw = model$zw.new - model$b.new * sqrt(model$w.new) - (n/2) * lambda0 * model$cw.new
  init.theta = rep(1, d)
  len = length(lambda_theta)
  measure <- matrix(0, ncol = len, nrow = nfolds)
  miss <- matrix(0, ncol = len, nrow = nfolds)

  for (f in 1:nfolds) {
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    # generate SS-ANOVA
    tr_n = length(trainID)
    te_n = length(testID)

    tr_G = G[trainID,]
    te_G = G[testID,]

    tr_RS = RiskSet(time[trainID], status[trainID])
    te_RS = RiskSet(time[testID], status[testID])
    for (k in 1:len) {
      # init.theta = as.vector(glmnet(Gw[trainID,], uw[trainID], family = "gaussian", lambda = lambda_theta[k])$beta)
      fit = gettheta.cd(init.theta, tr_G, time[trainID], status[trainID], model$b.new, (tr_n/2) * lambda0 * model$cw.new[trainID], lambda_theta[k], gamma, tr_RS)

      Gw = tr_G * sqrt(fit$w.new)
      XX = fit$z.new - tr_G %*% fit$theta.new - model$b.new
      num = t(XX) %*% diag(fit$w.new) %*% XX
      den = (1 - sum(diag( Gw %*% ginv( t(Gw) %*% Gw) %*% t(Gw) )) / tr_n)^2
      measure[f, k] <- as.vector(num / den / tr_n)

      Lik = PartialLik(time[trainID], status[trainID], tr_RS, tr_G, fit$theta.new, neg = FALSE)
      miss[f, k] = Lik
    }
  }
  print(measure)
  print(miss)
  measure[measure == -Inf | measure == Inf | is.nan(measure)] <- NA
  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22

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
  main = "Cox Family"
  max_min <- c(min(cvm - cvsd, na.rm = TRUE), max(cvm + cvsd, na.rm = TRUE))
  ylab = expression("GCV(" * lambda[theta] * ")")

  xrange = log(lambda_theta)
  plot(xrange, cvm, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab, ylim = max_min, type = 'n')
  arrows(xrange, cvm - cvsd, xrange, cvm + cvsd, angle = 90, code = 3, length = 0.1, col = 'gray')
  points(xrange, cvm, pch = 15, col = 'red')
  abline(v = xrange[id], col = 'darkgrey')
  # text(log(lambda_theta), par("usr")[4], labels = selm, pos = 1)
  if(one.std) abline(v = xrange[std.id], col = 'darkgrey')

  ##
  cvm <- apply(miss, 2, mean, na.rm = T)
  cvsd <- apply(miss, 2, sd, na.rm = T) / sqrt(nrow(miss)) + 1e-22
  max_min <- c(min(cvm - cvsd, na.rm = TRUE), max(cvm + cvsd, na.rm = TRUE))

  plot(log(lambda_theta), cvm, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = "partial likelihood", ylim = max_min, type = 'n')
  try(arrows(log(lambda_theta), cvm - cvsd, log(lambda_theta), cvm + cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(lambda_theta), cvm, pch = 15, col = 'red')
  abline(v = log(lambda_theta)[id], col = 'darkgrey', lty = 2)
  if(one.std) abline(v = log(lambda_theta)[std.id], col = 'darkgrey')

  if(algo == "CD"){
    # theta.new = .Call("Cnng", Gw, uw, n, d, init.theta, optlambda, gamma)
    # init.theta = as.vector(glmnet(Gw, uw, family = "gaussian", lambda = optlambda)$beta)
    fit = gettheta.cd(init.theta, G, time, status, model$b.new, (n/2) * lambda0 * model$cw.new, optlambda, gamma, model$RS)
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = fit$theta.new)
  }

  return(out)
}

# G = tr_G
# bhat = model$b.new
# const = (tr_n/2) * lambda0 * model$cw.new[trainID]
# lambda_theta = lambda_theta[1]
# gamma
# Risk = tr_RS
gettheta.cd = function(init.theta, G, time, status, bhat, const, lambda_theta, gamma, Risk){
  n = nrow(G)
  d = ncol(G)
  r = lambda_theta * gamma * n

  wz = calculate_wz_for_theta(init.theta, G, time, status, Risk)
  w = wz$weight
  z = wz$z

  uw = (z * sqrt(w)) - bhat * sqrt(w) - const
  Gw = G * sqrt(w)
  # uw = u * sqrt(w)
  theta = rep(1, d)
  theta.new <- temp <- init.theta
  for(iter in 1:20){
    for(j in 1:d){
      temp[j] = 2 * sum((uw - Gw[,-j] %*% theta[-j]) * Gw[,j])
      temp[j] = ifelse(temp[j] > 0 & r < abs(temp[j]), temp[j], 0)
      temp[j] = temp[j] / (sum(Gw[,j]^2) + n * lambda_theta * (1-gamma)) / 2

      loss = abs(theta - temp)
      conv1 = max(loss) < 1e-6
      # conv2 = sum(exp(G %*% theta.new) == Inf) > 0

      if(conv1) break
      theta[j] <- theta.new[j] <- temp[j]
    }
    if(conv1) break
  }
print(theta.new)
  # if(i == 1 & (conv1 | conv2)) theta.new = rep(0, d)

  return(list(z.new = z, w.new = w, theta.new = theta.new))
}


gettheta.cd = function(init.theta, G, time, status, bhat, const, lambda_theta, gamma, Risk){
  n = nrow(G)
  d = ncol(G)
  r = lambda_theta * gamma * n

  wz = calculate_wz_for_theta(init.theta, G, time, status, Risk)
  w = wz$weight
  z = wz$z

  uw = (z * sqrt(w)) - bhat * sqrt(w) - const
  Gw = G * sqrt(w)
  # uw = u * sqrt(w)
  theta = init.theta
  theta.new = rep(0, d)
  for(i in 1:20){
    for(j in 1:d){
      theta.new[j] = 2 * sum((uw - Gw[,-j] %*% theta[-j]) * Gw[,j])
      theta.new[j] = ifelse(theta.new[j] > 0 & r < abs(theta.new[j]), theta.new[j], 0)
      theta.new[j] = theta.new[j] / (sum(Gw[,j]^2) + n * lambda_theta * (1-gamma)) / 2

      loss = abs(theta - theta.new)
      conv = max(loss) < 1e-6

      if(conv) break
      theta[j] = theta.new[j]
    }
    if(conv) break
  }

  return(list(z.new = z, w.new = w, theta.new = theta.new))
}

calculate_wz_for_theta = function(init.theta, G, time, status, RS){
  n = length(time)
  weight = z = rep(0, n)

  for (k in 1:n) {
    Sum.exp.eta.Grad = Sum.exp.eta.Hess = 0
    id = which(RS[k,] > 0)
    eta = as.numeric(G[k,] %*% init.theta)
    exp.eta = exp(eta)
    for(r in id){
      Sum.exp.eta = sum(exp(G[RS[,r],] %*% init.theta))
      Sum.exp.eta.Grad = Sum.exp.eta.Grad + exp.eta / Sum.exp.eta # {j in R_i} exp(R_j c)
      Sum.exp.eta.Hess = Sum.exp.eta.Hess + ( exp.eta * Sum.exp.eta - exp.eta^2 ) / Sum.exp.eta^2
    }

    Grad.Term = status[k] - Sum.exp.eta.Grad
    weight[k] = Sum.exp.eta.Hess
    z[k] = eta + Grad.Term / weight[k]
  }

  return(list(z = z, weight = weight))
}


# calculate_GH_for_theta = function(theta, G, chat, time, status, lambda0, Risk){
#   n = length(time)
#   d = length(theta)
#   tie.size = as.numeric(table(time[status == 1]))
#   # if (min(eigen(R)$value) < 0)
#   #   R = R + 1e-08 * diag(nrow(R))
#   eta = G %*% theta
#   Hess.FullNumer.unScale = array(NA, dim = c(d, d, n))
#   for (i in 1:n) Hess.FullNumer.unScale[, , i] = G[i, ] %*% t(G[i, ])
#
#   Grad.Term1 = -t(G) %*% status/n
#   Grad.Term2 = matrix(NA, ncol = ncol(Risk), nrow = d)
#   Grad.Term3 = 2 * lambda0 * t(G) %*% chat
#   Grad.FullNumer = t(G) %*% diag(as.numeric(exp(eta)))
#   Grad.FullDenom = Hess.FullDenom = exp(eta)
#   Hess.FullNumer = Hess.FullNumer.unScale * array(rep(exp(eta), each = d^2),
#                                                   dim = c(d, d, n)
#   )
#   Hess.Term1 = Hess.Term2 = array(NA, dim = c(d, d, ncol(Risk)))
#   k = 1
#   tempSum.exp.eta = sum(exp(eta[Risk[, k]]), na.rm = TRUE)
#   temp.Gradient.numer = apply(Grad.FullNumer[, Risk[, k]], 1, sum, na.rm = TRUE)
#   temp.Hessian.numer = apply(Hess.FullNumer[, , Risk[, k]], c(1, 2), sum, na.rm = TRUE)
#   Grad.Term2[, k] = tie.size[k] * temp.Gradient.numer/tempSum.exp.eta
#   Hess.Term1[, , k] = temp.Hessian.numer/tempSum.exp.eta
#   Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
#   for (k in 2:ncol(Risk)) {
#     excludeID = Risk[, k - 1][!Risk[, k - 1] %in% Risk[, k]]
#     tempSum.exp.eta = tempSum.exp.eta - sum(exp(eta[excludeID]))
#     if (length(excludeID) > 1) {
#       temp.Gradient.numer = temp.Gradient.numer - apply(Grad.FullNumer[, excludeID], 1, sum)
#       temp.Hessian.numer = temp.Hessian.numer - apply(Hess.FullNumer[, , excludeID], c(1, 2), sum)
#     }
#     else {
#       temp.Gradient.numer = temp.Gradient.numer - Grad.FullNumer[, excludeID]
#       temp.Hessian.numer = temp.Hessian.numer - Hess.FullNumer[, , excludeID]
#     }
#     Grad.Term2[, k] = tie.size[k] * temp.Gradient.numer/tempSum.exp.eta
#     Hess.Term1[, , k] = temp.Hessian.numer/tempSum.exp.eta
#     Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
#   }
#   Grad.Term2 = apply(Grad.Term2, 1, sum)/n
#   Gradient = as.vector(Grad.Term1 + Grad.Term2 + Grad.Term3)
#   Hessian = apply(Hess.Term1, c(1, 2), sum)/n - apply(Hess.Term2, c(1, 2), sum)/n
#   return(list(Gradient = Gradient, Hessian = Hessian))
# }

