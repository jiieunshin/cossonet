RiskSet = function (time, status)
{ n = length(time)
  time = unlist(as.vector(time))
  status = unlist(as.vector(status))

  RiskSet = matrix(0, ncol = n, nrow = n)
  for (k in 1:n) {
    risk.id = which(time >= time[k])
    RiskSet[risk.id, k] = risk.id
  }
  return(RiskSet)
}

# time = unlist(y[,"time"])
# status = unlist(y["status"])
# mscale = rep(1, d)/wt^2
# nfolds = 5
# cand.lambda = exp(seq(log(2^{-30}), log(2^{4}), length.out = 40))
# lambda0 =  exp(seq(log(2^{-30}), log(2^{4}), length.out = 40))[1]
cv.getc = function(x, time, status, mscale, nfolds, lambda0, one.std, type, kparam, algo)
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

  c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox', lambda = cand.lambda[1], alpha = 0)$beta)
  c.init = c(scale(c.init))
  f.init = c(Rtheta %*% c.init)

  measure <- matrix(NA, ncol = length(cand.lambda), nrow = nfolds)
  for (f in 1:nfolds) {
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    # generate SS-ANOVA
    tr_n = length(trainID)
    te_n = length(testID)

    tr_RS = RiskSet(time[trainID], status[trainID])
    te_RS = RiskSet(time[testID], status[testID])

    tr_R = array(NA, c(tr_n, n, d))
    te_R = array(NA, c(te_n, n, d))

    for(j in 1:d){
      tr_R[, , j] = K$K[[j]][trainID, ]
      te_R[, , j] = K$K[[j]][testID, ]
    }

    tr_Rtheta <- wsGram(tr_R, mscale)
    te_Rtheta <- wsGram(te_R, mscale)

    for (k in 1:length(cand.lambda)) {
      # dyn.load("src/coxfuncs.dll")
      .Call("Cget_c", tr_Rtheta, Rtheta, n, tr_n, tr_RS, c.init, cand.lambda[k])
      if(algo == "CD"){
        c.new = getc.cd(tr_Rtheta, Rtheta, time[trainID], status[trainID], tr_RS, c.init, cand.lambda[k])
      }

      if(algo == "QP"){
        fit = getc.QP(tr_Rtheta, Rtheta, time[trainID], status[trainID], tr_RS, cand.lambda[k])

      }

      expRc = sapply(1:te_n, function(j) {
        id = unique(te_RS[,j])
        sum(exp(te_Rtheta[id,] %*% c.new))
      })

    PL = -c(t(status[testID]) %*% te_Rtheta %*% c.new) + t(status[testID]) %*% expRc + cand.lambda[k] * c(t(c.new) %*% Rtheta %*% c.new)
    measure[f, k] <- PL
    # pen1 = sum(diag(diag(status[testID]) %*% Rtheta[testID,] %*% ginv(sspline_fit$hessian) %*% t(Rtheta[testID,]) %*% diag(status[testID]))) / length(testID) / (length(testID)-1)
    # pen2 = sum(t(status[testID]) %*% Rtheta[testID,] %*% ginv(sspline_fit$hessian) %*% t(Rtheta[testID,]) %*% status[testID]) / length(testID)^2 / (length(testID)-1)
    # measure[f, k] <- PL + pen1 - pen2

    }
  }
  # print(measure)
  rm(tr_Rtheta)
  rm(te_Rtheta)
  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22

  ylab = "Partial likelihood"

  # optimal lambda1
  id = which.min(cvm)[1]
  optlambda = cand.lambda[id]

  # plotting error bar
  main = "Cox family"

  max_min <- c(min(cvm - cvsd), max(cvm + cvsd))

  plot(log(cand.lambda), cvm, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab, ylim = max_min, type = 'n')
  try(arrows(log(cand.lambda), cvm - cvsd, log(cand.lambda), cvm + cvsd, angle = 90, length = 0.01, col = 'gray'), silent = TRUE)
  points(log(cand.lambda), cvm, pch = 15, col = 'red')
  abline(v = log(cand.lambda)[id], col = 'darkgrey', lty = 2)
  # if(one.std) abline(v = log(cand.lambda)[std.id], col = 'darkgrey')

  if(algo == "CD"){
    c.new = getc.cd(Rtheta, Rtheta, time, status, RS, c.init, cand.lambda[k])
    f.new = c(Rtheta %*% c.new)
    out = list(IDmat = IDmat, measure = measure, R = R, RS = RS, f.new = f.new, c.new = c.new, optlambda = optlambda)
  }

  if(algo == "QP"){
    c.new = getc.QP(Rtheta, Rtheta, time, status, RS, c.init, cand.lambda[k])
    f.new = c(Rtheta %*% c.new)
    out = list(IDmat = IDmat, measure = measure, R = R, RS = RS, f.new = f.new, c.new = c.new, optlambda = optlambda)
  }
  rm(K)
  rm(Rtheta)

  return(out)
}

# R = tr_Rtheta
# time = time[trainID]
# status = status[trainID]
# Risk = tr_RS
# lambda0 = cand.lambda[1]

getc.cd = function(R1, R2, time, status, Risk, c.init, lambda0){
  m = nrow(R1)
  n = ncol(R1) # full basis

  # initialize
  c.old = c.init

  c.new = rep(0, length(c.old))
  for(iter in 1:80){

    for(k in 1:n){
      expRc = sapply(1:m, function(j) {
        id = unique(Risk[,j])
        sum(exp(R1[id,] %*% c.old))
      })
      print(expRc)
      RexpRc = sapply(1:m, function(j) {
        id = unique(Risk[,j])
        sum(R1[id,k] * exp(R1[id,] %*% c.old))
      })

      RRexpRc = sapply(1:m, function(j) {
        id = unique(Risk[,j])
        sum(R1[id,k]^2 * exp(R1[id,] %*% c.old))
      })

      grad = - sum(status * (R1[,k] - RexpRc / expRc)) + m * lambda0 * sum(c.old * R2[k,])
      hess = sum(status * (RRexpRc / expRc - (RexpRc * RexpRc) / (expRc * expRc) )) + m * lambda0 * R2[k,k]

      c.new[k] = c.old[k] - grad / hess

      loss = abs(c.old - c.new)
      conv = max(loss) < 1e-6

      if(conv) break
      c.old = c.new  # if not convergence

    } # end inner iteration
  } # end outer iteration

  if(iter == 1 & !conv) c.new = c.old

  return(c.new)
}

# model = sspline_cvfit
# mscale = wt
# lambda_theta = exp(seq(log(2^{-6}), log(2^{12}), length.out = 40))
# gamma = 0.8
cv.gettheta = function (model, x, time, status, mscale, lambda0, lambda_theta, gamma, nfolds, one.std, algo)
{
  n = dim(Kmat)[1]
  d = dim(Kmat)[3]
  IDmat = model$IDmat

  # solve theta
  G <- matrix(0, nrow(model$R[, ,1]), d)
  for (j in 1:d) {
    G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
  }

  theta_list = list()
  measure <- matrix(NA, ncol = length(lambda_theta), nrow = nfolds)
  l = 0
  for (f in 1:nfolds) {
    testID <- IDmat[!is.na(IDmat[, f]), f]
    trainID <- (1:n)[-testID]

    tr_G = G[trainID,]
    te_G = G[testID,]

    tr_n = length(trainID)
    te_n = length(testID)
    for (k in 1:length(lambda_theta)) {
      if(algo == "CD") {
        theta.fit = gettheta.cd(tr_G, G, model$chat, time[trainID], status[trainID], Risk[trainID,], model$optlambda, lambda_theta[k])
      }

      if(algo == "QP") {
        theta.fit = gettheta.QP(tr_G, G, model$chat, time[trainID], status[trainID], Risk[trainID,], model$optlambda, lambda_theta[k])
      }

      expG = sapply(1:nrow(Risk[testID,]), function(j) {exp(G[testID,] %*% theta.fit$thetahat)})
      expG = sum(log(rowSums(expG))) / n
      # print(theta.fit$thetahat)
      PL = -c(t(status[testID]) %*% G[testID, ] %*% theta.fit$thetahat) / length(testID) + expG +
        model$optlambda * c(t(model$chat[testID]) %*% G[testID, ] %*% theta.fit$thetahat) +
        lambda_theta[k] * gamma * sum(theta.fit$thetahat) + lambda_theta[k] * (1-gamma) * norm(theta.fit$thetahat, "2")
      measure[f, k] <- PL
      # pen1 = sum(diag(diag(status[testID]) %*% G[testID,] %*% ginv(theta.fit$hessian + diag(1e-8, d)) %*% t(G[testID,]) %*% diag(status[testID]))) / length(testID) / (length(testID)-1)
      # pen2 = sum(t(status[testID]) %*% G[testID,] %*% ginv(theta.fit$hessian + diag(1e-8, d)) %*% t(G[testID,]) %*% status[testID]) / length(testID)^2 / (length(testID)-1)
      # measure[f, k] <- PL + pen1 - pen2
      # print(measure[f, k])

    }
  }
  measure[is.nan(measure)] = 1e-30
  ylab = "Partial likelihood"

  cvm <- apply(measure, 2, mean, na.rm = T)
  cvsd <- apply(measure, 2, sd, na.rm = T) / sqrt(nrow(measure)) + 1e-22
  # selm = floor(apply(sel, 2, mean))

  id = which.min(cvm)[1]

  if(one.std){
    st1_err = cvm[id] + cvsd[id] # minimum cv err
    std.id = max(which(cvm[id:d] <= st1_err & cvm[id] <= cvm[id:d]))
    std.id = ifelse(std.id > id, std.id, id)
    optlambda = lambda_theta[std.id]
  } else{
    optlambda = lambda_theta[id]
  }

  # plotting error bar
  main = "Cox Family"
  max_min <- c(min(cvm - cvsd), max(cvm + cvsd))

  xrange = log(lambda_theta)
  plot(xrange, cvm, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = "generalized cross validation", ylim = max_min, type = 'n')
  arrows(xrange, cvm - cvsd, xrange, cvm + cvsd, angle = 90, code = 3, length = 0.1, col = 'gray')
  points(xrange, cvm, pch = 15, col = 'red')
  abline(v = xrange[id], col = 'darkgrey')
  # text(log(lambda_theta), par("usr")[4], labels = selm, pos = 1)
  if(one.std) abline(v = xrange[std.id], col = 'darkgrey', lty = 2)

  if(algo == "CD"){
    theta.new = gettheta.cd(G, G, model$chat, time, status, Risk, model$optlambda, optlambda)
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta_list = theta_list, thetahat = theta.fit$thetahat)
  }

  if(algo == "QP"){
    theta.new = gettheta.QP(G, G, model$chat, time, status, Risk, model$optlambda, optlambda)
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta_list = theta_list, thetahat = theta.fit$thetahat)
  }
  return(out)
}


# G1 = G1[trainID, ]
# G2
# time = time[trainID]
# status[trainID]
# Risk[trainID,]
# lambda_theta[k]

gettheta.cd = function(G1, G2, chat, time, status, Risk, lambda0, lambda_theta)
{
  n = nrow(G1)
  d = ncol(G1)

  theta.old = as.vector(coxph(Surv(time, status) ~ G1)$coefficients)
  theta.old = theta.old / sd(theta.old)
  theta.old[is.na(theta.old)] = 1e-6
  theta.old <- pmax(1e-6, theta.old)
  # print(theta.old)
  f.old = c(G1 %*% theta.old)

  theta.new = theta.old

  for(iter in 1:40){
    # for(j in 1:d){

    expG = sapply(1:nrow(Risk), function(j) {exp(t(G1[j,]) %*% theta.old)})
    GexpG_list = lapply(1:nrow(Risk), function(j) {G1[j,] * c(exp(t(G1[j,]) %*% theta.old))})
    GGexpG_list = lapply(1:nrow(Risk), function(j) {G1[j,] %*% t(G1[j,]) * c(exp(t(G1[j,]) %*% theta.old))})

    GexpG = matrix(0, d, 1)
    for(k in 1:nrow(Risk)) GexpG = GexpG + GexpG_list[[k]]

    GGexpG = matrix(0, d, d)
    for(k in 1:nrow(Risk)) GGexpG = GGexpG + GGexpG_list[[k]]

    grad = - t(G1) %*% status / n + GexpG / sum(expG) / n + lambda0 * t(G2) %*% chat
    hessian = (GGexpG / sum(expG) - (GexpG / sum(expG)) %*% t(GexpG / sum(expG)) ) / n

    upd = hessian %*% theta.old + grad + 2 * lambda_theta * (1-gamma) * theta.old + lambda_theta * gamma

    theta.new <- pmax(1e-6, theta.old - 0.01 * pmin(max(theta.old), upd))

    if(sum(theta.new) == 0 | sum(is.nan(theta.new)) > 0) theta.new <- theta.old
    if(sum(theta.new == Inf) >= 1) break
    if(max(abs(theta.old-theta.new)) <= 1e-4) break

    f.new = c(G1 %*% theta.new)
    theta.old = theta.new
  }

  if(sum(theta.new == Inf) >= 1) theta.new = theta.old
  # theta.new <- theta.new
  id = which(theta.new <= 1e-5)
  theta.new[id] <- 0
  # print(theta.new)
  out = list(iteration = iter, hessian = hessian, lambda_theta = lambda_theta, gamma = gamma, thetahat = theta.new)
  return(out)
}
