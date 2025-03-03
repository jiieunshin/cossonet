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

cv.getc.subset = function(K, time, status,  nbasis, basis.id, mscale, cand.lambda, type, nfold, kparam, one.std, show)
{
  cat("-- c-step -- \n")
  cat("proceeding... \n\n")

  d = K$numK
  n <- length(status)
  len = length(cand.lambda)

  Uv = array(NA, c(n, nbasis, d))
  for(j in 1:d){
    Uv[, , j] = K$K[[j]][, basis.id]
  }
  U <- combine_kernel(Uv, mscale)

  Qv = array(NA, c(nbasis, nbasis, d))
  for(j in 1:d){
    Qv[, , j] = K$K[[j]][basis.id, basis.id]
  }
  Q <- combine_kernel(Qv, mscale)

  # cross-validation
  fold = cvsplitID(n, nfold, status, family = "gaussian")
  measure <- matrix(NA, nfold, len)
  for(fid in 1:nfold){
    tr_id = as.vector(fold[, -fid])
    te_id = fold[, fid]

    tr_id = tr_id[!is.na(tr_id)]
    te_id = te_id[!is.na(te_id)]

    tr_n = length(tr_id)
    te_n = length(te_id)

    tr_U = array(NA, c(tr_n, nbasis, d))
    for(j in 1:d){
      tr_U[, , j] = K$K[[j]][tr_id, basis.id]
    }

    tr_U <- combine_kernel(tr_U, mscale)

    te_U = array(NA, c(te_n, nbasis, d))
    for(j in 1:d){
      te_U[, , j] = K$K[[j]][te_id, basis.id]
    }

    te_U <- combine_kernel(te_U, mscale)

    # initialize
    loop = 0
    EigQ = eigen(Q)
    while (min(eigen(Q)$values) < 0 & loop < 10) {
      loop = loop + 1
      Q = Q + 1e-08 * diag(nbasis)
      EigQ = eigen(Q)
    }
    if (loop == 10)
      EigQ$values[EigQ$values < 0] = 1e-08
    pseudoX = tr_U %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))

    for (k in 1:len){
      response <- survival::Surv(time = time[tr_id], event = status[tr_id])
      c.init = as.vector(glmnet(pseudoX, response, family = "cox", lambda = cand.lambda[k], alpha = 1, standardize = FALSE)$beta)
      eta = exp(tr_U %*% c.init)
      coxgrad_results <- coxgrad(eta, response, rep(1, tr_n), std.weights = FALSE, diag.hessian = TRUE)
      w <- - attributes(coxgrad_results)$diag_hessian
      z <- (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)

      zw = z * sqrt(w)
      Uw = tr_U * w
      sw = sqrt(w)

      fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, tr_n, nbasis, tr_n * cand.lambda[k], PACKAGE = "cossonet")
      c.new = fit$c.new

      testf = c(te_U %*% c.new)

        err = te_n * sum(w * (time[te_id] - testf)^2)
        inv.mat = ginv(t(te_U) %*% te_U + cand.lambda[k] * Q)
        df = sum(diag(te_U %*% inv.mat %*% t(te_U)))
        measure[fid, k] = err / (te_n - df)^2

      # if(cv == "ACV"){
      #   te_RS = RiskSet(time[te_id], status[te_id])
      #   tr_RS = RiskSet(time[tr_id], status[tr_id])
      #   test_GH = cosso::gradient.Hessian.C(fit$c.new, te_R, Q, time[te_id], status[te_id], mscale, cand.lambda[k], te_RS)
      #   train_GH = cosso::gradient.Hessian.C(fit$c.new, tr_R, Q, time[tr_id], status[tr_id], mscale, cand.lambda[k], tr_RS)
      #
      #   UHU = tr_U %*% My_solve(train_GH$H, t(tr_U))
      #   ACV_pen = sum(status[tr_id] == 1)/tr_n^2 * (sum(diag(UHU))/(tr_n - 1) - sum(UHU)/(tr_n^2 - tr_n))
      #   measure[fid, k] = PartialLik(time[tr_id], status[tr_id], tr_RS, tr_U %*% fit$c.new) + ACV_pen
      # }

    }
  }

  # smoothing parameter selection
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

  rm(tr_U)
  rm(te_U)

  RS = RiskSet(time, status)

  response = survival::Surv(time = time, event = status)
  c.init = as.vector(glmnet(U, response, family = "cox", lambda = optlambda, alpha = 1, standardize = FALSE)$beta)
  eta = exp(U %*% c.init)
  coxgrad_results <- coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
  w = - attributes(coxgrad_results)$diag_hessian
  z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)

  zw = z * sqrt(w)
  Uw = U * w
  sw = sqrt(w)

  fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, n, nbasis, n * optlambda, PACKAGE = "cossonet")

  out = list(measure = measure, Uv = Uv, Q = Q, RS = RS, f.new = c(U %*% fit$c.new),
             zw.new = zw, w.new = w, sw.new = sw, c.new = fit$c.new,
             optlambda = optlambda)

  rm(K)
  rm(Uv)
  rm(U)
  rm(Qv)
  rm(Q)
  rm(RS)
  return(out)
}

cv.gettheta.subset = function (model, K, time, status, nbasis, basis.id, mscale, lambda0, lambda_theta, gamma, nfold, one.std)
  {
  cat("-- theta-step -- \n")
  cat("proceeding... \n\n")

  n = length(time)
  d = length(mscale)

  G <- matrix(0, n, d)
  for (j in 1:d) {
    G[, j] = (model$Uv[, , j] %*% model$c.new) * (mscale[j]^(-2))
  }

  Gw <- matrix(0, n, d)
  for (j in 1:d) {
    Gw[, j] = ((model$Uv[, , j] * sqrt(model$w.new)) %*% model$c.new) * (mscale[j]^(-2))
  }

  uw = model$zw.new - model$sw.new

  h = rep(0, d)
  for (j in 1:d) {
    h[j] = n * ((t(model$c.new) %*% model$Uv[basis.id, , j]) %*% model$c.new) * lambda0
  }

  # cross-validation
  init.theta = rep(1, d)
  len = length(lambda_theta)
  measure = matrix(NA, nfold, len)
  fold = cvsplitID(n, nfold, time, family = "gaussian")

  for(fid in 1:nfold){
    tr_id = as.vector(fold[, -fid])
    te_id = fold[, fid]

    tr_id = tr_id[!is.na(tr_id)]
    te_id = te_id[!is.na(te_id)]

    tr_n = length(tr_id)
    te_n = length(te_id)

    te_Uv = array(NA, c(te_n, nbasis, d))
    for(j in 1:d){
      te_Uv[, , j] = K$K[[j]][te_id, basis.id]
    }

    for (k in 1:len) {
      response <- survival::Surv(time = time[tr_id], event = status[tr_id])
      eta = exp(G[tr_id,] %*% init.theta)
      coxgrad_results = coxgrad(eta, response, rep(1, tr_n), std.weights = FALSE, diag.hessian = TRUE)
      w = - attributes(coxgrad_results)$diag_hessian
      z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0) + lambda0 * G[tr_id,] %*% t(G[basis.id, ]) %*% model$c.new

      theta.new = .Call("wls_theta_step", Gw[tr_id,], uw[tr_id], h/2, tr_n, d, init.theta, tr_n * lambda_theta[k] * gamma / 2, tr_n * lambda_theta[k] * (1-gamma), PACKAGE = "cossonet")
      theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)

      te_U = wsGram(te_Uv, theta.adj/mscale^2)
      fhat = c(te_U %*% model$c.new)

      err = te_n * sum(w * (time[te_id] - fhat)^2)
      inv.mat = ginv(t(te_U) %*% te_U + lambda_theta[k] * model$Q)
      df = sum(diag(te_U %*% inv.mat %*% t(te_U)))
      measure[fid, k] = err / (te_n - df)^2

      # if(cv == "ACV") {
      #   ACV = cosso::PartialLik(time[tr_id], status[tr_id], RiskSet(time[tr_id], status[tr_id]), fhat) + model$ACV_pen
      #   measure[fid, k] = ACV
      # }

    }
  }

  rm(te_Uv)
  rm(te_U)

  # smoothing parameter selection
  measure_mean = colMeans(measure, na.rm = T)
  measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(5)

  sel_id = which(!is.nan(measure_se) & measure_se != Inf)
  measure_mean = measure_mean[sel_id]
  measure_se = measure_se[sel_id]
  lambda_theta = lambda_theta[sel_id]

  min_id = which.min(measure_mean)


  if(one.std){
    cand_ids = which((measure_mean >= measure_mean[min_id]) &
                       (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
    cand_ids = cand_ids[cand_ids >= min_id]
    std_id = max(cand_ids)
    optlambda = lambda_theta[std_id]
  } else{
    optlambda = lambda_theta[min_id]
  }

  ylab = expression("GCV(" * lambda[theta] * ")")

  plot(log(lambda_theta), measure_mean, main = "Cox family", xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab,
       ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
  arrows(x0 = log(lambda_theta), y0 = measure_mean - measure_se,
         x1 = log(lambda_theta), y1 = measure_mean + measure_se,
         angle = 90, code = 3, length = 0.1, col = "darkgray")
  abline(v = log(optlambda), lty = 2, col = "darkgray")

  response = survival::Surv(time = time, event = status)
  eta = exp(G %*% init.theta)
  coxgrad_results <- coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
  w = - attributes(coxgrad_results)$diag_hessian
  z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0) + lambda0 * G %*% t(G[basis.id, ]) %*% model$c.new

  theta.new = .Call("wls_theta_step", Gw, uw, h/2, n, d, init.theta, n * optlambda * gamma / 2, n * optlambda * (1-gamma), PACKAGE = "cossonet")
  theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)

  out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.adj)

  rm(G)
  rm(Gw)

  return(out)
}

soft_threshold = function(a, b){
  return(ifelse(a > 0 & b < abs(a), a - b, 0))
}
