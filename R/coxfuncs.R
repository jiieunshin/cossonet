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

# mscale = wt
# cand.lambda = lambda0
cv.getc = function(K, time, status, mscale, cand.lambda, type, kparam, algo, show)
{
  d = K$numK
  n <- length(time)
  len = length(cand.lambda)

  R = array(NA, c(n, n, d))
  for(j in 1:d){
    R[, , j] = K$K[[j]]
  }

  Rtheta <- combine_kernel(R, mscale)

  f.init = rep(0.5, n)
  RS = RiskSet(time, status)

  measure <- rep(0, length(cand.lambda))
  gcv_list  <- rep(0, length(cand.lambda))
  for (k in 1:length(cand.lambda)){
    if(algo == "CD"){
      EigRtheta = eigen(Rtheta)
      if (min(EigRtheta$value) < 0) {
        Rtheta = Rtheta + max(1e-07, 1.5 * abs(min(EigRtheta$value))) * diag(nrow(Rtheta))
        EigRtheta = eigen(Rtheta)
      }
      pseudoX = Rtheta %*% EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values))
      ssCox.en = glmnet(pseudoX, cbind(time = time, status = status),
                        family = "cox", lambda = cand.lambda, alpha = 0,
                        standardize = FALSE)
      c.init = as.numeric(EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values)) %*% ssCox.en$beta[, 1])

      # f.init = c(Rtheta %*% c.init)
      fit = getc.cd(R, Rtheta, mscale, f.init, c.init, time, status, cand.lambda[k], RS)

      # Rw = Rtheta * fit$c.new
      # XX = fit$zw.new - Rw %*% fit$cw.new - fit$b.new * sqrt(fit$w.new)
      # num = t(XX) %*% XX + 1
      # S = Rw %*% ginv(t(Rw) %*% Rw) %*% t(Rw)
      # den = (1 - sum(diag(S)) / n)^2 + 1
      # measure[k] = as.vector( num / den / n )

      measure[k] = fit$GCV

      # gcv_list[k] = fit$GCV
    }

    if(algo == "QP"){
      # c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
      #                           lambda = cand.lambda[k], alpha = 0)$beta)

      EigRtheta = eigen(Rtheta)
      if (min(EigRtheta$value) < 0) {
        Rtheta = Rtheta + max(1e-07, 1.5 * abs(min(EigRtheta$value))) * diag(nrow(Rtheta))
        EigRtheta = eigen(Rtheta)
      }
      pseudoX = Rtheta %*% EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values))
      ssCox.en = glmnet(pseudoX, cbind(time = time, status = status),
                        family = "cox", lambda = cand.lambda, alpha = 0,
                        standardize = FALSE)
      c.init = as.numeric(EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values)) %*% ssCox.en$beta[, 1])

      fit = getc.QP(R, Rtheta, c.init, time, status, mscale, cand.lambda[k], RS)

      measure[k] = fit$GCV
    }
  }

  id = which.min(measure)[1]
  optlambda = cand.lambda[id]

  # optimal lambda1
  if(show) plot(log(cand.lambda), measure, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = "partial likelihood", ylim = range(measure), pch = 15, col = 'red')
  # if(show) plot(log(cand.lambda), gcv_list, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = "partial likelihood", ylim = range(measure), pch = 15, col = 'red')

  if(algo == "CD"){
    EigRtheta = eigen(Rtheta)
    if (min(EigRtheta$value) < 0) {
      Rtheta = Rtheta + max(1e-07, 1.5 * abs(min(EigRtheta$value))) * diag(nrow(Rtheta))
      EigRtheta = eigen(Rtheta)
    }
    pseudoX = Rtheta %*% EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values))
    ssCox.en = glmnet(pseudoX, cbind(time = time, status = status),
                      family = "cox", lambda = optlambda, alpha = 0,
                      standardize = FALSE)
    c.init = as.numeric(EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values)) %*% ssCox.en$beta[, 1])

    # f.init = c(Rtheta %*% c.init)
    fit = getc.cd(R, Rtheta, mscale, f.init, c.init, time, status, optlambda, RS)
    out = list(measure = measure, R = R, f.new = c(Rtheta %*% fit$c.new),
               w.new = fit$w.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)
    }

  if(algo == "QP"){

    EigRtheta = eigen(Rtheta)
    if (min(EigRtheta$value) < 0) {
      Rtheta = Rtheta + max(1e-07, 1.5 * abs(min(EigRtheta$value))) * diag(nrow(Rtheta))
      EigRtheta = eigen(Rtheta)
    }
    pseudoX = Rtheta %*% EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values))
    ssCox.en = glmnet(pseudoX, cbind(time = time, status = status),
                      family = "cox", lambda = optlambda, alpha = 0,
                      standardize = FALSE)
    c.init = as.numeric(EigRtheta$vectors %*% diag(sqrt(1/EigRtheta$values)) %*% ssCox.en$beta[, 1])

    # c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
    #                           lambda = optlambda, alpha = 0, standardize = FALSE)$beta)
    fit = getc.QP(R, Rtheta, c.init, time, status, mscale, optlambda, RS)
    out = list(measure = measure, R = R, f.new = c(Rtheta %*% fit$c.new), c.new = fit$c.new, optlambda = optlambda, conv = TRUE)
  }

  rm(K)
  rm(Rtheta)

  return(out)
}

getc.cd = function(R, Rtheta, mscale, f, c.init, time, status, lambda0, Risk)
{
  n = ncol(Rtheta)
  # wz = calculate_wz_for_c(c.init, Rtheta, time, status, Risk)
  # w = wz$weight
  # z = wz$z

  # return(list(zw.new = zw, w.new = w, sw.new = sw, b.new = b.new, c.new = c.new, cw.new = cw.new))

  c.old = c.init
  c.new = rep(0, n)
  # while (loop < 15 & iter.diff > 1e-4) {

    for(i in 1:10){ # outer iteration
      GH = try(calculate_GH_for_C(c.old, R, R, time, status, mscale, lambda0, Risk), silent = TRUE)
      err = class(GH) == "try-error"
      if(err) break
      Hess = GH$Hessian
      Grad = GH$Gradient
      # 2 * n * lambda0 * Rtheta2
      W = ginv(Hess)
      z = (Hess %*% c.old - Grad) / lambda0
      for(j in 1:n){
        V1 = t(z - Rtheta[ ,-j] %*% c.old[-j]) %*% (t(W) %*% Rtheta[, j])
        V2 = (c.old[-j] %*% Rtheta[-j, j]) / lambda0
        V3 = t(Rtheta[, j]) %*% (t(W) %*% Rtheta[, j])
        V4 = Rtheta[j, j] / lambda0

        c.new[j] = (V1 - V2) / (V3 + V4)
        loss = abs(c.old - c.new) / abs(c.old)
        conv1 = min(loss[loss > 0]) < 1e-4
        conv2 = abs(c.old[j] - c.new[j]) > 5
        # cat("i = ", i, "j = ", j, "loss =", max(loss),  "\n")
        if(conv1 | conv2) break
        c.old[j] = c.new[j]  # if not convergence
      }
      if(conv1 | conv2 | err) break
    }

    if(i == 1 & (conv1 | conv2 | err)) c.new = c.init

  # zw = z * sqrt(w)
  # Rw = Rtheta * w
  # cw = c.init
  # cw.new = temp = c.init / sqrt(w)
  # sw = sqrt(w)
  # fit = .Call("cox_c_step", zw, Rw, cw, sw, n, lambda0, PACKAGE = "cdcosso")
  #
  # b.new = fit$b.new
  # c.new = fit$c.new
  # cw.new = fit$cw.new

  z = (Hess %*% c.new - Grad) / lambda0
  loglik = t(z - Rtheta %*% c.new) %*% W %*% (z - Rtheta %*% c.new)
  den = (1 - sum(diag(Rtheta %*% ginv(Rtheta + Hess/lambda0))) / n)^2
  GCV = as.numeric(loglik / den / n)
# print(i)

  # UHU = Rtheta %*% My_solve(GH$H, t(Rtheta))
  # ACV = PartialLik(time, status, Risk, Rtheta %*% c.new) + sum(status == 1)/n^2 * (sum(diag(UHU))/(n - 1) - sum(UHU)/(n^2 - n))
  return(list(z.new = z, w.new = W, c.new = c.new, GCV = GCV))
  # return(list(z.new = z, zw.new = zw, w.new = w, c.new = c.new, b.new = b.new, cw.new = cw.new, GCV = GCV))
}

calculate_GH_for_C = function (initC, Gramat1, Gramat2, time, status, mscale, lambda0, riskset, Hess.FullNumer.unScale)
{
  n = length(time)
  tie.size = as.numeric(table(time[status == 1]))
  Rtheta1 = wsGram(Gramat1, mscale)
  Rtheta2 = wsGram(Gramat2, mscale)
  if (min(eigen(Rtheta2)$value) < 0)
    Rtheta2 = Rtheta2 + 1e-08 * diag(nrow(Rtheta2))
  eta = Rtheta1 %*% initC
  if (missing(Hess.FullNumer.unScale)) {
    Hess.FullNumer.unScale = array(NA, dim = c(length(initC), length(initC), n))
    for (i in 1:n) Hess.FullNumer.unScale[, , i] = Rtheta1[i, ] %*% t(Rtheta1[i, ])
  }
  Grad.Term1 = -t(Rtheta1) %*% status
  Grad.Term2 = matrix(NA, ncol = ncol(riskset), nrow = length(initC))
  # Grad.Term3 = 2 * n * lambda0 * Rtheta2 %*% initC
  Grad.FullNumer = t(Rtheta1) %*% diag(as.numeric(exp(eta)))
  Grad.FullDenom = Hess.FullDenom = exp(eta)
  Hess.FullNumer = Hess.FullNumer.unScale * array(rep(exp(eta),
                                                      each = length(initC)^2),
                                                  dim = c(length(initC), length(initC), n)
                                                  )
  Hess.Term1 = Hess.Term2 = array(NA, dim = c(length(initC), length(initC), ncol(riskset)))
  k = 1
  tempSum.exp.eta = sum(exp(eta[riskset[, k]]), na.rm = TRUE)
  temp.Gradient.numer = apply(Grad.FullNumer[, riskset[, k]], 1, sum, na.rm = TRUE)
  temp.Hessian.numer = apply(Hess.FullNumer[, , riskset[, k]], c(1, 2), sum, na.rm = TRUE)
  Grad.Term2[, k] = tie.size[k] * temp.Gradient.numer/tempSum.exp.eta
  Hess.Term1[, , k] = temp.Hessian.numer/tempSum.exp.eta
  Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
  for (k in 2:ncol(riskset)) {
    excludeID = riskset[, k - 1][!riskset[, k - 1] %in% riskset[, k]]
    tempSum.exp.eta = tempSum.exp.eta - sum(exp(eta[excludeID]))
    if (length(excludeID) > 1) {
      temp.Gradient.numer = temp.Gradient.numer - apply(Grad.FullNumer[, excludeID], 1, sum)
      temp.Hessian.numer = temp.Hessian.numer - apply(Hess.FullNumer[, , excludeID], c(1, 2), sum)
    }
    else {
      temp.Gradient.numer = temp.Gradient.numer - Grad.FullNumer[, excludeID]
      temp.Hessian.numer = temp.Hessian.numer - Hess.FullNumer[, , excludeID]
    }
    Grad.Term2[, k] = tie.size[k] * temp.Gradient.numer/tempSum.exp.eta
    Hess.Term1[, , k] = temp.Hessian.numer/tempSum.exp.eta
    Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
  }
  Grad.Term2 = apply(Grad.Term2, 1, sum)
  Gradient = Grad.Term1 / n + Grad.Term2 / n
  Hessian = apply(Hess.Term1, c(1, 2), sum) / n - apply(Hess.Term2, c(1, 2), sum) / n
  return(list(Gradient = Gradient, Hessian = Hessian))
}

getc.QP = function (R, Rtheta, c.init, time, status, mscale, lambda0, RS)
{
  n = length(time)
  p = length(mscale)

  GH = cosso::gradient.Hessian.C(c.init, R, R, time, status, mscale, lambda0, RS)
  c.new = as.numeric(cosso::My_solve(GH$H, GH$H %*% c.init - GH$G))
  UHU = Rtheta %*% cosso::My_solve(GH$H, t(Rtheta))

  W = ginv(GH$H)
  z = (GH$H %*% c.new - GH$G) / lambda0

  loglik = t(z - Rtheta %*% c.new) %*% W %*% (z - Rtheta %*% c.new)
  den = (1 - sum(diag(Rtheta %*% ginv(Rtheta + GH$H/lambda0))) / n)^2
  GCV = as.numeric(loglik / den / n)
  return(list(c.new = c.new, GCV = GCV))
}

calculate_wz_for_c = function(c.init, R, time, status, RS){
  n = length(time)
  Grad.Term = weight = z = rep(0, n)

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

    Grad.Term[k] = status[k] - Sum.exp.eta.Grad
    weight[k] = Sum.exp.eta.Hess
    z[k] = eta - ifelse(weight[k] != 0, - Grad.Term[k]/weight[k], 0)
  }

  return(list(z = z, gradient = Grad.Term, weight = weight))
}

# model = getc_cvfit
# lambda0 = getc_cvfit$optlambda
# mscale = wt
cv.gettheta = function (model, x, time, status, mscale, lambda0, lambda_theta, gamma, type, kparam, algo){
  n = length(time)
  d = length(mscale)
  IDmat = model$IDmat

  RS = RiskSet(time, status)

  # solve theta
  G <- matrix(0, nrow(model$R[, ,1]), d)
  for (j in 1:d) {
    G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
  }

  if(algo == "QP") lambda_theta = exp(seq(log(1e-4), log(40), length.out = length(lambda_theta)))
  len = length(lambda_theta)

  measure <- rep(0, len)
  save_theta <- list()
  for (k in 1:len) {
    if(algo == "CD"){
      # init.theta = rep(1, d)

      # Gw = G * sqrt(model$w.new)
      # uw = model$zw.new - model$b.new * sqrt(model$w.new) - (n/2) * lambda0 * model$cw.new
      # theta.new = .Call("cox_theta_step", Gw, uw, n, d, rep(1, d), lambda_theta[k], gamma)
      # save_theta[[k]] <- theta.new
      #
      # XX = model$zw.new - Gw %*% theta.new
      # num = t(XX) %*% XX + 1
      # den = (1 - sum(diag( Gw %*% ginv( t(Gw) %*% Gw) %*% t(Gw) )) / n)^2 + 1

      # fit = gettheta.cd(rep(1, d), model$f.new, G, time, status, model$b.new * sqrt(model$w.new), model$c.new,
      #                   (1/2) * lambda0 * model$cw.new,
      #                   lambda0, lambda_theta[k], gamma, RS)
      fit = gettheta.cd(rep(1, d), model$f.new, G, time, status, 0, model$c.new,
                        0, lambda0, lambda_theta[k], gamma, RS)

      save_theta[[k]] <- fit$theta.new

      # num = as.vector(cosso::PartialLik(time, status, RS,  G %*% fit$theta.new)) + 1
      # den = (1 - sum(fit$theta.new != 0))^2 + 1
      # measure[k] <- num / den /n

      measure[k] <- fit$ACV
    }

    if(algo == "QP"){
      fit = gettheta.QP(rep(1, d), model$c.new, G, time, status, lambda0, lambda_theta[k], RS)
      save_theta[[k]] <- fit$theta.new

      # measure[k] <- cosso::PartialLik(time, status, RS,  G %*% fit$theta.new) + sum(status == 1)/n^2 * (sum(diag(fit$UHU))/(n - 1) - sum(fit$UHU)/(n^2 - n))
    }
  }
  # print(save_theta)
  # print(measure)
  id = which.min(measure)[1]
  optlambda = lambda_theta[id]

  # plotting error bar
  xrange = log(lambda_theta)
  plot(xrange, measure, main = "Cox family", xlab = expression("Log(" * lambda[theta] * ")"), ylab = "partial likelihood", ylim = range(measure), pch = 15, col = 'red')

  if(algo == "CD"){
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = save_theta[[id]])
    # fit = gettheta.cd(init.theta, G, time, status, model$b.new, (n/2) * lambda0 * model$cw.new, optlambda, gamma, RS)
    # out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = fit$theta.new)
  }
  if(algo == "QP"){
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = save_theta[[id]])
    # fit = gettheta.QP(init.theta, model$c.new, G, time, status, lambda0, optlambda, RS)
    # out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = fit$theta.new)
  }

  return(out)
}

gettheta.cd = function(init.theta, f.init, G, time, status, bhat, chat, const, lambda0, lambda_theta, gamma, Risk){
  n = nrow(G)
  d = ncol(G)
  r = lambda_theta * gamma

  # wz = calculate_wz_for_theta(theta.old, G, time, status, Risk)
  # w = wz$weight
  # z = wz$z

  Hess.FullNumer.unScale = array(NA, dim = c(length(init.theta), length(init.theta), n))
  for (i in 1:n) Hess.FullNumer.unScale[, , i] = G[i, ] %*% t(G[i, ])

  theta.old = init.theta
  theta.new = rep(0, d)
  for(i in 1:20){
    GH = cosso::gradient.Hessian.Theta(theta.old, chat, G, G, lambda0, 0, time, status, Risk, Hess.FullNumer.unScale)
    Dmat = GH$Hessian / 2
    dvec = Dmat %*% theta.old - GH$Gradient
    for(j in 1:d){
      # if(j == 1){
      #   L = 0
      #   U = Hess[1, 2:d] %*% theta.old[2:d]
      # } else if(j == d){
      #   L = Hess[d, 1:(d-1)] %*% theta.old[1:(d-1)]
      #   U = 0
      # } else{
      #   L = Hess[j, 1:(j-1)] %*% theta.old[1:(j-1)]
      #   U = Hess[j, (j+1):d] %*% theta.old[(j+1):d]
      # }
      theta.new[j] = dvec[j] - Dmat[j, -j] %*% theta.old[-j]
        # Hess[j, -j] %*% theta.old[-j]
      theta.new[j] = soft_threshold(theta.new[j], r)
      theta.new[j] = theta.new[j] / (Dmat[j, j] + lambda_theta * (1-gamma))

      loss = abs(theta.old - theta.new) / abs(theta.old + 1)
      # cat("i = ", i, "j =", j, "theta.new[j] =", theta.new[j], "loss =", max(loss), "\n")
      conv = max(loss) < 1e-6

      if(conv) break
      theta.old[j] = theta.new[j]
    }
    if(conv) break
  }

  if(i == 1 & !conv) theta.new = rep(0, d)

  # Hess.FullNumer.unScale = array(NA, dim = c(length(init.theta), length(init.theta), n))
  # for (i in 1:n) Hess.FullNumer.unScale[, , i] = G[i, ] %*% t(G[i, ])
  #
  # GH = gradient.Hessian.Theta(init.theta, chat, G, G, lambda0, time, status, Risk, Hess.FullNumer.unScale)
  # Hess = GH$Hessian
  # Grad = GH$Gradient
  # # 2 * n * lambda0 * Rtheta2
  # w = 1/diag(Hess)
  # z = f.init - ifelse(w != 0, - Grad/w, 0)
  #
  # uw = c(z * sqrt(w)) - bhat - const
  # Gw = G * sqrt(w)
  #
  # theta.new = .Call("cox_theta_step", Gw, uw, n, d, init.theta, lambda_theta, gamma)
  # theta.new = ifelse(theta.new <= 1e-6, 0, theta.new)

  UHU = G %*% My_solve(GH$H, t(G))
  ACV = cosso::PartialLik(time, status, Risk,  G %*% theta.new) + sum(status == 1)/n^2 * (sum(diag(UHU))/(n - 1) - sum(UHU)/(n^2 - n))

  return(list(theta.new = theta.new, ACV = ACV))

  # return(list(Gw = Gw, zw.new = z * sqrt(w), uw.new = uw, w.new = w, theta.new = theta.new))
}

soft_threshold = function(a, b){
 return(sign(a) * ifelse(a > 0 & b < abs(a), a - b, 0))
}

gradient.Hessian.Theta = function (initTheta, initC, G1, G2, lambda0, time, status, riskset, Hess.FullNumer.unScale)
{
  n = length(time)
  p = length(initTheta)
  tie.size = as.numeric(table(time[status == 1]))
  eta = G1 %*% initTheta
  Grad.Term1 = -t(G1) %*% status/n
  Grad.Term2 = matrix(NA, ncol = ncol(riskset), nrow = p)
  Grad.Term3 = lambda0 * t(G2) %*% initC / 2
  Grad.FullNumer = t(G1) %*% diag(as.numeric(exp(eta)))
  Grad.FullDenom = Hess.FullDenom = exp(eta)
  Hess.FullNumer = Hess.FullNumer.unScale * array(rep(exp(eta), each = p^2), dim = c(p, p, n))
  Hess.Term1 = Hess.Term2 = array(NA, dim = c(p, p, ncol(riskset)))
  k = 1
  tempSum.exp.eta = sum(exp(eta[riskset[, k]]), na.rm = TRUE)
  tempGradient.numer = apply(Grad.FullNumer[, riskset[, k]], 1, sum, na.rm = TRUE)
  tempHessian.numer = apply(Hess.FullNumer[, , riskset[, k]], c(1, 2), sum, na.rm = TRUE)
  Grad.Term2[, k] = tie.size[k] * tempGradient.numer/tempSum.exp.eta
  Hess.Term1[, , k] = tempHessian.numer/tempSum.exp.eta
  Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
  for (k in 2:ncol(riskset)) {
    excludeID = riskset[, k - 1][!riskset[, k - 1] %in% riskset[, k]]
    tempSum.exp.eta = tempSum.exp.eta - sum(exp(eta[excludeID]))
    if (length(excludeID) > 1) {
      tempGradient.numer = tempGradient.numer - apply(Grad.FullNumer[, excludeID], 1, sum)
      tempHessian.numer = tempHessian.numer - apply(Hess.FullNumer[, , excludeID], c(1, 2), sum)
    }
    else {
      tempGradient.numer = tempGradient.numer - Grad.FullNumer[, excludeID]
      tempHessian.numer = tempHessian.numer - Hess.FullNumer[, , excludeID]
    }
    Grad.Term2[, k] = tie.size[k] * tempGradient.numer/tempSum.exp.eta
    Hess.Term1[, , k] = tempHessian.numer/tempSum.exp.eta
    Hess.Term2[, , k] = 1/tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
  }
  Grad.Term2 = apply(Grad.Term2, 1, sum)/n
  Gradient = Grad.Term1 + Grad.Term2 + Grad.Term3
  Hessian = apply(Hess.Term1, c(1, 2), sum)/n - apply(Hess.Term2, c(1, 2), sum)/n
  return(list(Gradient = Gradient, Hessian = Hessian))
}

calculate_wz_for_theta = function(init.theta, G, time, status, RS){
  n = length(time)
  Grad.Term = weight = z = rep(0, n)

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

    Grad.Term[k] = status[k] - Sum.exp.eta.Grad
    weight[k] = Sum.exp.eta.Hess
    z[k] = eta - ifelse(weight[k] != 0, - Grad.Term[k]/weight[k], 0)
  }

  return(list(z = z, gradient = Grad.Term, weight = weight))
}

gettheta.QP = function(init.theta, c.hat, G, time, status, lambda0, lambda_theta, Risk){
  n = nrow(G)
  p = ncol(G)
  Hess.FullNumer.unScale = array(NA, dim = c(length(init.theta),
                                             length(init.theta),
                                             n)
                                 )
  for (i in 1:n) Hess.FullNumer.unScale[, , i] = G[i, ] %*% t(G[i, ])
  loop = 0
  iter.diff = Inf
  old.Theta = init.theta
  while (loop < 15 & iter.diff > 1e-04) {
    loop = loop + 1
    GH = cosso::gradient.Hessian.Theta(old.Theta, c.hat, G, G,
                                       lambda0, lambda_theta, time, status, Risk, Hess.FullNumer.unScale)

    if(min(eigen(GH$H)$value) < 0)
      GH$H = GH$H + max(1e-07, 1.5 * abs(min(eigen(GH$H)$value))) * diag(length(old.Theta))

    dvec = -(GH$G - GH$H %*% old.Theta)
    Amat = t(rbind(diag(p), rep(-1, p)))
    bvec = c(rep(0, p), -lambda_theta)
    new.Theta = cosso::My_solve.QP(GH$H, dvec, Amat, bvec)
    new.Theta[new.Theta < 1e-07] = 0
    iter.diff = mean(abs(new.Theta - old.Theta))
    old.Theta = new.Theta
  }

  UHU = G %*% My_solve(GH$H, t(G))

  return(list(theta.new = new.Theta, G = GH$G, H = GH$H, UHU = UHU))
}
