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
cv.getc.subset = function(K, time, status, nbasis, basis.id, mscale, cand.lambda, type, kparam, one.std, show)
{
  d = K$numK
  n <- length(status)
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

  fold = cvsplitID(n, 5, status, family = "gaussian")
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

    # initialize
    loop = 0
    EigRtheta2 = eigen(Rtheta2)
    while (min(eigen(Rtheta2)$values) < 0 & loop < 10) {
      loop = loop + 1
      Rtheta2 = Rtheta2 + 1e-08 * diag(nbasis)
      EigRtheta2 = eigen(Rtheta2)
    }
    if (loop == 10)
      EigRtheta2$values[EigRtheta2$values < 0] = 1e-08
    pseudoX = tr_Rtheta %*% EigRtheta2$vectors %*% diag(sqrt(1/EigRtheta2$values))

    for (k in 1:len){
      response <- survival::Surv(time = time[tr_id], event = status[tr_id])
      c.init = as.vector(glmnet(tr_Rtheta, response, family = "cox", lambda = cand.lambda[k], alpha = 1, standardize = FALSE)$beta)
      eta = exp(tr_Rtheta %*% c.init)
      coxgrad_results <- coxgrad(eta, response, rep(1, tr_n), std.weights = FALSE, diag.hessian = TRUE)
      w <- - attributes(coxgrad_results)$diag_hessian
      z <- (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
      fit <- .Call("cox_c_step", c.init, tr_Rtheta, Rtheta2, as.integer(tr_n), as.integer(nbasis), z, w, cand.lambda[k])

      # fit = getc.cd(tr_R, R2, tr_Rtheta, Rtheta2, mscale, c.init, time[tr_id], status[tr_id], cand.lambda[k], tr_RS)

      # calculate ACV for test data
      te_RS = RiskSet(time[te_id], status[te_id])

      test_GH = cosso::gradient.Hessian.C(fit$c.new, te_R, R2, time[te_id], status[te_id], mscale, cand.lambda[k], te_RS)

#       test_GH = .Call("gradient_Hessian_C", fit$c.new, as.integer(tr_n), as.integer(nbasis), as.integer(ncol(te_RS)), exp(te_Rtheta %*% fit$c.new),
#                       te_Rtheta, Rtheta2, time[te_id], as.integer(status[te_id]), mscale, cand.lambda[k], as.integer(te_RS),
#                       as.integer(table(time[te_id][status[te_id] == 1])),
#                       PACKAGE = "cdcosso")

      UHU = te_Rtheta %*% My_solve(test_GH$H, t(te_Rtheta))
      ACV_pen = sum(status[te_id] == 1)/te_n^2 * (sum(diag(UHU))/(te_n - 1) - sum(UHU)/(te_n^2 - te_n))

      measure[f, k] = PartialLik(time[te_id], status[te_id], te_RS, te_Rtheta %*% fit$c.new) + ACV_pen
    }
  }

  # optimal lambda1
  measure_mean = colMeans(measure, na.rm = T)
  measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(5)

  sel_id = which(!is.nan(measure_se) & measure_se != Inf)
  measure_mean = measure_mean[sel_id]
  measure_se = measure_se[sel_id]
  cand.lambda = cand.lambda[sel_id]

  min_id = which.min(measure_mean)

  if(one.std){
    cand_ids = which((measure_mean >= measure_mean[min_id]) &
                       (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
    cand_ids = cand_ids[cand_ids >= min_id]
    std_id = max(cand_ids)
    optlambda = cand.lambda[std_id]
  } else{
    optlambda = cand.lambda[min_id]
  }


  ylab = expression("GCV(" * lambda[0] * ")")

  if(show){
    plot(log(cand.lambda), measure_mean, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
         ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
    arrows(x0 = log(cand.lambda), y0 = measure_mean - measure_se,
           x1 = log(cand.lambda), y1 = measure_mean + measure_se,
           angle = 90, code = 3, length = 0.1, col = "darkgray")
    abline(v = log(optlambda), lty = 2, col = "darkgray")
  }

  rm(tr_R)
  rm(te_R)
  rm(tr_Rtheta)
  rm(te_Rtheta)
  rm(test_GH)
  rm(te_RS)
  # rm(tr_RS)

  response <- survival::Surv(time = time, event = status)
  c.init = as.vector(glmnet(Rtheta, response, family = "cox", lambda = optlambda, alpha = 1, standardize = FALSE)$beta)
  eta = exp(Rtheta %*% c.init)
  coxgrad_results <- coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
  w <- - attributes(coxgrad_results)$diag_hessian
  z <- (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
  fit <- .Call("cox_c_step", c.init, Rtheta, Rtheta2, as.integer(n), as.integer(nbasis), z, w, optlambda)


  # c.init = as.vector(glmnet(pseudoX, cbind(time, status), family = "cox", lambda = optlambda, alpha = 1, standardize = FALSE)$beta)
  RS = RiskSet(time, status)
  # fit = getc.cd(R, R2, Rtheta, Rtheta2, mscale, c.init, time, status, optlambda, RS)

  GH = cosso::gradient.Hessian.C(fit$c.new, R, R2, time, status, mscale, optlambda, RS)

  # GH =  .Call("gradient_Hessian_C", fit$c.new, as.integer(n), as.integer(nbasis), as.integer(ncol(RS)), exp(Rtheta %*% fit$c.new),
  #             R, R2, time, as.integer(status), mscale, optlambda, as.integer(RS), as.integer(table(time[status == 1])),
  #             PACKAGE = "cdcosso")

  UHU = Rtheta %*% My_solve(GH$H, t(Rtheta))
  ACV_pen = sum(status == 1)/n^2 * (sum(diag(UHU))/(n - 1) - sum(UHU)/(n^2 - n))

  out = list(measure = measure, R = R, RS = RS, f.new = c(Rtheta %*% fit$c.new), w.new = fit$w.new,
             c.new = fit$c.new, ACV_pen = ACV_pen, optlambda = optlambda)

  rm(K)
  rm(Rtheta)
  rm(Rtheta2)
  rm(GH)
  rm(UHU)
  return(out)
}

# R = tr_R
# Rtheta = tr_Rtheta
# time = time[tr_id]
# status = status[tr_id]
# lambda0 = cand.lambda[k]
# Risk = tr_RS
# getc.cd(tr_R, R2, tr_Rtheta, Rtheta2, mscale, c.init, time[tr_id], status[tr_id], cand.lambda[k], tr_RS)

getc.cd = function(R, R2, Rtheta, Rtheta2, mscale, c.init, time, status, lambda0, Risk)
{
  n = nrow(Rtheta)
  m = ncol(Rtheta)
  c.old = c.init
  c.new = rep(0, m)
  # GH = .Call("gradient_Hessian_C", c.old, as.integer(n), as.integer(m), as.integer(ncol(Risk)), as.numeric(exp(Rtheta %*% c.old)),
  #            as.numeric(Rtheta), Rtheta2, as.numeric(time), as.integer(status), mscale, lambda0, as.integer(Risk),
  #            as.integer(table(time[status == 1])),
  #            PACKAGE = "cdcosso")


  # while (loop < 15 & iter.diff > 1e-4) {
  for(i in 1:15){ # outer iteration
    # 2 * n * lambda0 * Rtheta2
    GH = cosso::gradient.Hessian.C(c.old, R, R2, time, status, mscale, lambda0, Risk)

    err = (class(GH) == "try-error") | sum(is.nan(GH$Gradient)) > 0
    if(err) break

    Hess = GH$Hessian - 2 * lambda0 * Rtheta2
    Grad = GH$Gradient - 2 * lambda0 * Rtheta2 %*% c.old

    W = ginv(Hess)
    z = (Hess %*% c.old - Grad) / lambda0

    for(j in 1:m){
      WR = colSums(W * Rtheta2[, j])
      V1 = sum(WR * (z - Rtheta2[ ,-j] %*% c.old[-j]))
      V2 = as.vector((Rtheta2[-j, j] %*% c.old[-j]) / lambda0)
      V3 = sum(WR * Rtheta2[ ,j])
      V4 = Rtheta2[j, j] / lambda0

      c.new[j] = (V1 - V2) / (V3 + V4)
      loss = abs(c.old - c.new)
      conv1 = min(loss[loss > 0]) < 1e-20
      conv2 = abs(c.old[j] - c.new[j]) > 5
      conv3 = sum(exp(Rtheta %*% c.new) == Inf) > 0
      # cat("i = ", i, "j = ", j, "loss =", max(loss),  "\n")
      if(conv1 | conv2 | conv3) break
      c.old[j] = c.new[j]  # if not convergence
    }
    if(conv1 | conv2 | conv3) break
  }

  if(i == 1 & (conv1 | conv2 | conv3)) c.new = c.init

  return(list(z.new = z, w.new = W, c.new = c.new))
  # return(list(z.new = z, zw.new = zw, w.new = w, c.new = c.new, b.new = b.new, cw.new = cw.new, GCV = GCV))
}


# model = getc_cvfit
# lambda0 = getc_cvfit$optlambda
# mscale = wt
cv.gettheta.subset = function (model, K, time, status, nbasis, basis.id, mscale, lambda0, lambda_theta, gamma){
  n = length(time)
  d = length(mscale)

  # RS = RiskSet(time, status)

  # solve theta
  G <- matrix(0, n, d)
  for (j in 1:d) {
    G[, j] = (model$R[, , j] %*% model$c.new) * (mscale[j]^(-2))
  }


  init.theta = rep(1, d)
  len = length(lambda_theta)
  measure = matrix(NA, 5, len)
  fold = cvsplitID(n, 5, time, family = "gaussian")

  for(f in 1:5){
    tr_id = as.vector(fold[, -f])
    te_id = fold[, f]

    tr_id = tr_id[!is.na(tr_id)]
    te_id = te_id[!is.na(te_id)]

    tr_n = length(tr_id)
    te_n = length(te_id)

    for (k in 1:len) {
      fit = gettheta.cd(rep(1, d), model$f.new[tr_id], G[tr_id, ], G[basis.id, ], time[tr_id], status[tr_id], model$c.new,
                        lambda0, lambda_theta[k], gamma, RiskSet(time[tr_id], status[tr_id]))

      # save_theta[[k]] = fit$theta.new

      ACV = cosso::PartialLik(time[te_id], status[te_id], RiskSet(time[te_id], status[te_id]), G[te_id, ] %*% fit$theta.new) + model$ACV_pen
      measure[f, k] = ACV

    }
  }
  print(measure)
  measure_mean = colMeans(measure, na.rm = T)
  measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(5)

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


  plot(log(lambda_theta), measure_mean, main = "Cox family", xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab,
       ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
  arrows(x0 = log(lambda_theta), y0 = measure_mean - measure_se,
         x1 = log(lambda_theta), y1 = measure_mean + measure_se,
         angle = 90, code = 3, length = 0.1, col = "darkgray")
  abline(v = log(lambda_theta)[std_id], lty = 2, col = "darkgray")

  fit = gettheta.cd(rep(1, d), model$f.new, G, G[basis.id, ], time, status, model$c.new,
                    lambda0, optlambda, gamma, RiskSet(time, status))

  theta.adj = ifelse(fit$theta.new <= 1e-6, 0, fit$theta.new)

  out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.adj)

  return(out)
}


# G1 = G[tr_id, ]
# G2 = G[basis.id,]
# time = time[tr_id]
# status = status[tr_id]
# init.theta = rep(1, d)
# f.init = model$f.new[tr_id]
# chat = model$c.new
# ACV_pen = model$ACV_pen
# lambda_theta = lambda_theta[k]
# Risk = RiskSet(time, status)
gettheta.cd = function(init.theta, f.init, G1, G2, time, status, chat, lambda0, lambda_theta, gamma, Risk){
  n = nrow(G1)
  d = ncol(G1)
  r = lambda_theta * gamma

  # wz = calculate_wz_for_theta(theta.old, G, time, status, Risk)
  # w = wz$weight
  # z = wz$z

  Hess.FullNumer.unScale = array(NA, dim = c(length(init.theta), length(init.theta), n))
  for (i in 1:n) Hess.FullNumer.unScale[, , i] = G1[i, ] %*% t(G1[i, ])

  theta.old = init.theta
  theta.new = rep(0, d)
  conv2 = conv3 = TRUE

  for(i in 1:20){
    GH = GH.theta(theta.old, chat, G1, G2, lambda0, time, status, Risk, Hess.FullNumer.unScale)

    loss = rep(1, d)
    err = sum(is.nan(GH$Gradient)) > 0
    if (err) break
    Dmat = GH$H / 2
    dvec = - (GH$H %*% theta.old - GH$Gradient)
    for(j in 1:d){
      if(j == 1){
        L = 0
        U = Dmat[1, 2:d] %*% theta.old[2:d]
      } else if(j == d){
        L = Dmat[d, 1:(d-1)] %*% theta.old[1:(d-1)]
        U = 0
      } else{
        L = Dmat[j, 1:(j-1)] %*% theta.old[1:(j-1)]
        U = Dmat[j, (j+1):d] %*% theta.old[(j+1):d]
      }

      theta.new[j] = soft_threshold(-dvec[j] - L + U, r)
      # L + U
      # Dmat[j, -j] %*% theta.old[-j]
      D_diag = ifelse(Dmat[j, j] <= 0, 0, Dmat[j, j])
      theta.new[j] = theta.new[j] / (D_diag + lambda_theta * (1-gamma))

      # loss = abs(theta.old - theta.new)
      # conv = max(loss) < 1e-12
      loss[j] = abs(theta.old[j] - theta.new[j])
      # print(theta.new)
      conv2 = sum(loss == 0) == d | is.infinite(theta.new[j]) | is.na(theta.new[j])
      # conv3 = max(loss) > 5

      # cat("i = ", i, "j =", j, "theta.new[j] =", theta.new[j], "loss =", max(loss), "\n")
      if(conv2){
        conv = TRUE
      } else{
        conv = max(loss[loss > 0]) < 1e-18
      }

      if(conv) break
      theta.old[j] = theta.new[j]
    }
    if(conv) break
  }
  print(i)
  # print(theta.new)
  if(i == 1 & (conv2)) theta.old = rep(0, d)

  return(list(theta.new = theta.old))

  # return(list(Gw = Gw, zw.new = z * sqrt(w), uw.new = uw, w.new = w, theta.new = theta.new))
}

soft_threshold = function(a, b){
  return(ifelse(a > 0 & b < abs(a), a - b, 0))
}

# initTheta = init.theta
# initC = chat
# G1 = G[tr_id,]
# G2 = G[basis.id,]
# riskset = RS
GH.theta = function (initTheta, initC, G1, G2, lambda0, time, status, riskset, Hess.FullNumer.unScale)
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
    } else {
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
