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

  for (k in 1:length(cand.lambda)){
    if(algo == "CD"){
      c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                                lambda = cand.lambda[k], alpha = 0)$beta)
      fit = getc.cd(Rtheta, f.init, c.init, time, status, cand.lambda[k], RS)

      Rw = Rtheta * fit$w.new
      XX = fit$zw.new - Rw %*% fit$cw.new - fit$b.new * sqrt(fit$w.new)
      num = t(XX) %*% XX + 1

      S = Rw %*% ginv(t(Rw) %*% Rw) %*% t(Rw)
      den = (1 - sum(diag(S)) / n)^2 + 1
      measure[k] <- as.vector( num / den / n )

      # W = outer(fit$gradient, fit$gradient)
      # UHU = Rtheta %*% W %*% t(Rtheta)

      # measure[k] <- cosso::PartialLik(time, status, RS, Rtheta %*% fit$c.new)
      # + sum(status == 1)/n^2 * (sum(diag(UHU))/(n - 1) - sum(UHU)/(n^2 - n))
    }

    if(algo == "QP"){
      c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                                lambda = cand.lambda[k], alpha = 0)$beta)
      fit = getc.QP(R, Rtheta, c.init, time, status, mscale, cand.lambda[k], RS)

      # measure[k] <- cosso::PartialLik(time, status, RS, Rtheta %*% fit$c.new) + sum(status == 1)/n^2 * (sum(diag(fit$UHU))/(n - 1) - sum(fit$UHU)/(n^2 - n))

      HH =  fit$H - 2 * cand.lambda[k] * Rtheta
      HHH = ginv(HH/cand.lambda[k] + Rtheta)
      GG = fit$G - 2 * cand.lambda[k] * Rtheta %*% fit$c.new

      z = (HHH %*% fit$c.new - GG) / cand.lambda[k]
      num = t(z - Rtheta %*% fit$c.new) %*% ginv(HH) %*% (z - Rtheta %*% fit$c.new)
      S = Rtheta %*% ginv(Rtheta + HH/cand.lambda[k])
      den = (1 - sum(diag(S)) / n)^2 + 1
      measure[k] <- as.vector( num / den / n )
    }
  }

  id = which.min(measure)[1]
  optlambda = cand.lambda[id]

  # optimal lambda1
  if(show) plot(log(cand.lambda), measure, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = "partial likelihood", ylim = range(measure), pch = 15, col = 'red')

  if(algo == "CD"){
    c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                              lambda = optlambda, alpha = 0, standardize = FALSE)$beta)
    fit = getc.cd(Rtheta, f.init, c.init, time, status, optlambda, RS)
    out = list(measure = measure, R = R, f.new = c(Rtheta %*% fit$c.new) + fit$b.new, zw.new = fit$zw.new, w.new = fit$w.new,
               b.new = fit$b.new, cw.new = fit$cw.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)
    }

  if(algo == "QP"){
    c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                              lambda = optlambda, alpha = 0, standardize = FALSE)$beta)
    fit = getc.QP(R, Rtheta, c.init, time, status, mscale, optlambda, RS)
    out = list(measure = measure, R = R, f.new = Rtheta %*% fit$c.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)
  }

  rm(K)
  rm(Rtheta)

  return(out)
}


getc.cd = function(Rtheta, f, c.init, time, status, lambda0, Risk)
{
  n = ncol(Rtheta)
  # wz = calculate_wz_for_c(c.init, Rtheta, time, status, Risk)
  # w = wz$weight
  # z = wz$z
  b = 0

  y = cbind(time = time, status = status)
  coxgrad_results = coxgrad(f, y, rep(1, length(f)), std.weights = FALSE, diag.hessian = TRUE)
  w = -attributes(coxgrad_results)$diag_hessian
  z = (f - b) - ifelse(w != 0, -coxgrad_results/w, 0)

  zw = z * sqrt(w)
  Rw = Rtheta * w
  cw = c.init
  cw.new = temp = c.init / sqrt(w)
  sw = sqrt(w)
  fit = .Call("c_step", zw, Rw, cw, sw, n, lambda0, PACKAGE = "cdcosso")

  b.new = fit$b.new
  c.new = fit$c.new
  cw.new = fit$cw.new

  return(list(zw.new = zw, w.new = w, sw.new = sw, b.new = b.new, c.new = c.new, cw.new = cw.new))
}

getc.QP = function (R, Rtheta, c.init, time, status, mscale, lambda0, RS)
{
  n = length(time)
  p = length(mscale)

  GH = cosso::gradient.Hessian.C(c.init, R, R, time, status, mscale, lambda0, RS)
  c.new = as.numeric(cosso::My_solve(GH$H, GH$H %*% c.init - GH$G))
  UHU = Rtheta %*% cosso::My_solve(GH$H, t(Rtheta))
  return(list(c.new = c.new, G = GH$G, H = GH$H, UHU = UHU))
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
    z[k] = eta + (Grad.Term[k] + 0.1) / (weight[k] + 0.1)
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
      init.theta = rep(1, d)

      Gw = G * sqrt(model$w.new)
      uw = model$zw.new - model$b.new * sqrt(model$w.new) - (n/2) * lambda0 * model$cw.new
      theta.new = .Call("theta_step", Gw, uw, n, d, init.theta, lambda_theta[k], gamma)
      save_theta[[k]] <- theta.new

      XX = model$zw.new - Gw %*% theta.new
      num = t(XX) %*% XX + 1
      den = (1 - sum(diag( Gw %*% ginv( t(Gw) %*% Gw) %*% t(Gw) )) / n)^2 + 1


      # fit = gettheta.cd(init.theta, model$f.new, G, time, status, model$b.new, (n/2) * lambda0 * model$cw.new, lambda_theta[k], gamma, RS)
      # save_theta[[k]] <- fit$theta.new
      #
      # XX = fit$zw.new - fit$Gw %*% fit$theta.new
      # num = t(XX) %*% XX + 1
      # den = (1 - sum(diag( fit$Gw %*% ginv( t(fit$Gw) %*% fit$Gw) %*% t(fit$Gw) )) / n)^2 + 1

      measure[k] <- as.vector(num / den / n)

      # measure[k] <- cosso::PartialLik(time, status, RS, G %*% theta.adj) / (1 - sum(fit$theta.new != 0) / n)^2 / n
    }

    if(algo == "QP"){
      init.theta = rep(1, d)
      fit = gettheta.QP(init.theta, model$c.new, G, time, status, lambda0, lambda_theta[k], RS)
      save_theta[[k]] <- fit$theta.new

      measure[k] <- as.vector(cosso::PartialLik(time, status, RS,  G %*% fit$theta.new) / (1 - sum(fit$theta.new != 0))^2 /n)

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

gettheta.cd = function(init.theta, f.init, G, time, status, bhat, const, lambda_theta, gamma, Risk){
  n = nrow(G)
  d = ncol(G)
  r = lambda_theta * gamma * n

  # wz = calculate_wz_for_theta(init.theta, G, time, status, Risk)
  # w = wz$weight
  # z = wz$z
  # f.init = rep(0.5, n)
  y = cbind(time = time, status = status)
  coxgrad_results = coxgrad(f.init, y, rep(1, nrow(G)), std.weights = FALSE, diag.hessian = TRUE)
  w = - attributes(coxgrad_results)$diag_hessian
  z = f.init - ifelse(w != 0, - coxgrad_results/w, 0)

  uw = (z * sqrt(w)) - bhat * sqrt(w) - const
  Gw = G * sqrt(w)

  theta.new = .Call("theta_step", Gw, uw, n, d, init.theta, lambda_theta, gamma)
  # theta.new = ifelse(theta.new <= 1e-6, 0, theta.new)
  return(list(Gw = Gw, zw.new = z * sqrt(w), w.new = w, theta.new = theta.new))
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
    z[k] = eta + (Grad.Term[k] + 0.1) / (weight[k] + 0.1)
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
