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

partial_liklihood = function (time, status, RS, K, a, neg = FALSE) {
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
cv.getc = function(x, time, status, mscale, cand.lambda, one.std, type, kparam, algo)
{
  n <- length(time)
  # IDmat <- cvsplitID(n, nfolds)

  K = make_anovaKernel(x, x, type = type, kparam)
  d = K$numK
  R = array(NA, c(n, n, d))

  for(j in 1:d){
    R[, , j] = K$K[[j]]
  }

  Rtheta <- wsGram(R, mscale)

  RS = RiskSet(time, status)

  measure <- miss <- rep(0, length(cand.lambda))

  for (k in 1:length(cand.lambda)){
    if(algo == "CD"){
      c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                                lambda = cand.lambda[k], alpha = 0, standardize = FALSE)$beta)
      # zw = z * sqrt(w)
      # Rw = Rtheta * w
      # cw = c.init / sqrt(w)
      # sw = sqrt(w)
      fit = getc.cd(Rtheta, c.init, time, status, cand.lambda[k], RS)
      # fit = .Call("Csspline", tr_Rtheta, Rtheta, n, n, RS, c.init, cand.lambda[k])
      Lik = partial_liklihood(time, status, RS, Rtheta, fit$c.new, neg = FALSE)

      Rw = Rtheta * fit$w.new
      XX = fit$zw.new - Rw %*% fit$cw.new - fit$b.new * fit$w.new
      num = t(XX) %*% XX
      den = (1 - sum(diag(Rtheta %*% ginv(Rtheta + diag(fit$w.new)/cand.lambda[k]))) / n)^2

      # measure[k] <- as.vector(num / den / n)
      measure[k] <- -Lik
      # miss[k] = -Lik
    }

    if(algo == "QP"){
      c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                                lambda = cand.lambda[k], alpha = 0, standardize = FALSE)$beta)
      fit = getc.QP(R, Rtheta, c.init, time, status, mscale, cand.lambda[k], RS)
      measure[k] <- cosso::PartialLik(time, status, RS, fit)
    }

  }

  id = which.min(measure)[1]
  optlambda = cand.lambda[id]

  # optimal lambda1
  plot(log(cand.lambda), measure, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = "partial likelihood", ylim = range(measure), pch = 15, col = 'red')

  plot(log(cand.lambda), miss, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = "miss", ylim = range(miss), pch = 15, col = 'red')


  if(algo == "CD"){
    c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                              lambda = optlambda, alpha = 0, standardize = FALSE)$beta)
    fit = getc.cd(Rtheta, c.init, time, status, optlambda, RS)
    out = list(measure = measure, R = R, zw.new = fit$zw.new, w.new = fit$w.new,
               b.new = fit$b.new, cw.new = fit$cw.new, c.new = fit$c.new, optlambda = optlambda, conv = TRUE)
    }

  if(algo == "QP"){
    c.init = as.vector(glmnet(Rtheta, cbind(time = time, status = status), family = 'cox',
                              lambda = optlambda, alpha = 0, standardize = FALSE)$beta)
    fit = getc.QP(R, Rtheta, c.init, time, status, mscale, optlambda, RS)

    out = list(measure = measure, R = R, c.new = fit, optlambda = optlambda, conv = TRUE)
  }

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

getc.QP = function (R, Rtheta, c.init, time, status, mscale, lambda0, RS)
{
  n = length(time)
  p = length(mscale)

  GH = cosso::gradient.Hessian.C(c.init, R, R, time, status, mscale, lambda0, RS)
  c.new = as.numeric(cosso::My_solve(GH$H, GH$H %*% c.init - GH$G))
  # UHU = Rtheta %*% My_solve(GH$H, t(Rtheta))
  # print(cw.new)
  # out = list(cw.new = cw.new)
  return(c.new)
}

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
cv.gettheta = function (model, x, time, status, mscale, lambda0, lambda_theta, gamma, one.std, type, kparam, algo){
  n = length(time)
  d = length(mscale)
  IDmat = model$IDmat

  RS = RiskSet(time, status)

  # solve theta
  G <- matrix(0, nrow(model$R[, ,1]), d)
  for (j in 1:d) {
    G[, j] = model$R[, , j] %*% model$c.new * (mscale[j]^(-2))
  }

  init.theta = rep(1, d)

  if(algo == "QP") lambda_theta = exp(seq(log(0.2), log(80), length.out = length(lambda_theta)))
  len = length(lambda_theta)

  measure <- miss <- rep(0, len)

  for (k in 1:len) {
    if(algo == "CD"){
      fit = gettheta.cd(init.theta, G, time, status, model$b.new, (n/2) * lambda0 * model$cw.new, lambda_theta[k], gamma, RS)

      Gw = G * sqrt(fit$w.new)
      XX = fit$z.new - G %*% fit$theta.new - model$b.new
      num = t(XX) %*% diag(fit$w.new) %*% XX
      den = (1 - sum(diag( Gw %*% ginv( t(Gw) %*% Gw) %*% t(Gw) )) / n)^2
      # measure[k] <- as.vector(num / den / n)
      measure[k] <- -Lik
      Lik = partial_liklihood(time, status, RS, G, fit$theta.new, neg = FALSE)
      # miss[k] = -Lik
    }

    if(algo == "QP"){
      fit = gettheta.QP(init.theta, model$c.new, G, time, status, lambda0, lambda_theta[k], RS)
      measure[k] <- cosso::PartialLik(time, status, RS, G %*% fit)
    }
  }
  id = which.min(measure)[1]
  optlambda = lambda_theta[id]

  # plotting error bar
  xrange = log(lambda_theta)
  plot(xrange, measure, main = "Cox family", xlab = expression("Log(" * lambda[theta] * ")"), ylab = "partial likelihood", ylim = range(measure), pch = 15, col = 'red')

  plot(xrange, miss, main = "Cox family", xlab = expression("Log(" * lambda[theta] * ")"), ylab = "miss", ylim = range(miss), pch = 15, col = 'red')

  if(algo == "CD"){
    # theta.new = .Call("Cnng", Gw, uw, n, d, init.theta, optlambda, gamma)
    # init.theta = as.vector(glmnet(Gw, uw, family = "gaussian", lambda = optlambda)$beta)
    fit = gettheta.cd(init.theta, G, time, status, model$b.new, (n/2) * lambda0 * model$cw.new, optlambda, gamma, RS)
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = fit$theta.new)
  }

  if(algo == "QP"){
    fit = gettheta.QP(init.theta, model$c.new, G, time, status, lambda0, optlambda, RS)
    out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = fit)
  }

  return(out)
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
      conv = max(loss) < 1e-6 | is.na(theta.new[j])

      if(conv) break
      theta[j] = theta.new[j]
    }
    if(conv) break
  }

  if(i == 1 & !conv) theta = rep(0, d)

  return(list(z.new = z, w.new = w, theta.new = theta))
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
      GH$H = GH$H + max(1e-07, 1.5 * abs(min(eigen(GH$H)$value))) * diag(length(init.theta))

    dvec = -(GH$G - GH$H %*% old.Theta)
    Amat = t(rbind(diag(p), rep(-1, p)))
    bvec = c(rep(0, p), -lambda_theta)
    new.Theta = cosso::My_solve.QP(GH$H, dvec, Amat, bvec)
    new.Theta[new.Theta < 1e-07] = 0
    iter.diff = mean(abs(new.Theta - old.Theta))
    old.Theta = new.Theta
  }
  return(new.Theta)
}
