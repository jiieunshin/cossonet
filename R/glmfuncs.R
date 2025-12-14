cv.sspline.subset = function (K, y, nbasis, basis.id, mscale, cand.lambda, obj, type, cv, nfold, kparam, one.std, show)
{
  message("-- c-step -- \n")
  message("proceeding... \n")
  d = K$numK
  n <- length(y)
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

  EigQ = eigen(Q)
  loop = 0
  while (min(EigQ$values) < 0 & loop < 10) {
    loop = loop + 1
    Q = Q + 1e-08 * diag(nbasis)
    EigQ = eigen(Q)
  }
  if (loop == 10)
    EigQ$values[EigQ$values < 0] = 1e-08

  # cv = "GCV" -> use only training data
  if(cv == "GCV"){
    pseudoX = U %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))
    measure = numeric(NA, len)
    for (k in 1:len){
      c.init = as.vector(glmnet(pseudoX, y, family = obj$family, lambda = cand.lambda[k], alpha = 1, standardize = FALSE)$beta)
      
      ff = U %*% c.init
      mu = obj$linkinv(ff)
      w = as.vector(obj$variance(mu))
      z = ff + (y - mu) / w
      
      zw = z * sqrt(w)
      Uw = U * w
      sw = sqrt(w)
      
      fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, n, nbasis, n * cand.lambda[k], PACKAGE = "cossonet")
      b.new = fit$b.new
      c.new = fit$c.new
      
      testf = c(b.new + U %*% c.new)
      testmu = obj$linkinv(testf)
      testw = as.vector(obj$variance(testmu))
      
      err = n * sum(testw * (y - testf)^2)
      inv.mat = ginv(t(U) %*% U + cand.lambda[k] * Q)
      df = sum(diag(U %*% inv.mat %*% t(U)))
      measure[k] = err / (n - df)^2
    }
      rm(inv.mat)
      
      # smoothing parameter selection
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
      min_id = which.min(measure_mean)
      optlambda = cand.lambda[min_id]
      
      if(show){
        plot(log(cand.lambda), measure, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
             ylim = range(measure), pch = 15, col = 'red', type = "b")
        abline(v = log(optlambda), lty = 2, col = "darkgray")
      }
  }

  # cv = "mse" -> nfold
  if(cv == "mse"){
    
    fold = cvsplitID(n, nfold, y, family = obj$family)
    measure = matrix(NA, nfold, len)
    for(fid in 1:nfold){
      tr_id = as.vector(fold[, -fid])
      te_id = fold[, fid]
  
      tr_id = tr_id[!is.na(tr_id)]
      te_id = te_id[!is.na(te_id)]
  
      tr_n = length(tr_id)
      te_n = length(te_id)
  
      tr_Uv = array(NA, c(tr_n, nbasis, d))
      for(j in 1:d){
        tr_Uv[, , j] = K$K[[j]][tr_id, basis.id]
      }
      tr_U <- combine_kernel(tr_Uv, mscale)
  
      te_Uv = array(NA, c(te_n, nbasis, d))
      for(j in 1:d){
        te_Uv[, , j] = K$K[[j]][te_id, basis.id]
      }
      te_U <- combine_kernel(te_Uv, mscale)
  
      pseudoX = tr_U %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))
  
      for (k in 1:len){
        c.init = as.vector(glmnet(pseudoX, y[tr_id], family = obj$family, lambda = cand.lambda[k], alpha = 1, standardize = FALSE)$beta)
  
        ff = tr_U %*% c.init
        mu = obj$linkinv(ff)
        w = as.vector(obj$variance(mu))
        z = ff + (y[tr_id] - mu) / w
  
        zw = z * sqrt(w)
        Uw = tr_U * w
        sw = sqrt(w)
  
        fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, tr_n, nbasis, tr_n * cand.lambda[k], PACKAGE = "cossonet")
        b.new = fit$b.new
        c.new = fit$c.new
  
        testf = c(b.new + te_U %*% c.new)
        testmu = obj$linkinv(testf)
        testw = as.vector(obj$variance(testmu))
  
        if(obj$family == "gaussian") measure[fid, k] <- mean((testf - y[te_id])^2)
        if(obj$family == "binomial") measure[fid, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
        if(obj$family == "poisson") measure[fid, k] <- mean((y[te_id] - testf)^2)
      }
      
      # smoothing parameter selection
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
      
      measure_mean = colMeans(measure, na.rm = T)
      measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(5)
      
      sel_id = which(!is.nan(measure_se) & measure_se != Inf)
      measure_mean = measure_mean[sel_id]
      measure_se = measure_se[sel_id]
      cand.lambda = cand.lambda[sel_id]
      
      min_id = which.min(measure)
      
      if(one.std){
        cand_ids = which((measure_mean >= measure_mean[min_id]) &
                           (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
        cand_ids = cand_ids[cand_ids >= min_id]
        std_id = max(cand_ids)
        optlambda = cand.lambda[std_id]
      } else{
        optlambda = cand.lambda[min_id]
      }
      
      if(show){
        plot(log(cand.lambda), measure_mean, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
             ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
        arrows(x0 = log(cand.lambda), y0 = measure_mean - measure_se,
               x1 = log(cand.lambda), y1 = measure_mean + measure_se,
               angle = 90, code = 3, length = 0.1, col = "darkgray")
        abline(v = log(optlambda), lty = 2, col = "darkgray")
      }
      
      # if(cv == "KL"){
      #   true_mu = obj$linkinv(f[te_id])
      #   measure[fid, k] <- KL(testf, true_mu, obj)
      # }
      
      rm(tr_U)
      rm(te_U)
    }
  }

  pseudoX = U %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))
  c.init = as.vector(glmnet(pseudoX, y, family = obj$family, lambda = optlambda, alpha = 1, standardize = FALSE)$beta)

  ff = U %*% c.init
  mu = obj$linkinv(ff)
  w = as.vector(obj$variance(mu))
  z = ff + (y - mu) / w

  zw = z * sqrt(w)
  Uw = U * w
  sw = sqrt(w)

  fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, n, nbasis, n * optlambda, PACKAGE = "cossonet")
  b.new = fit$b.new
  c.new = fit$c.new

  f.new = c(b.new + U %*% c.new)
  mu.new = obj$linkinv(f.new)
  w.new = obj$variance(mu.new)
  z.new = f.new + (y - mu.new) / w.new

  if(obj$family == "gaussian") m = mean((f.new - y)^2)
  if(obj$family == "binomial") m <- mean(y != ifelse(mu.new < 0.5, 0, 1))
  if(obj$family == "poisson") m <- mean((y - f.new)^2)

  message("mse:", round(m, 4), "\n")

  out = list(measure = measure, Uv = Uv, Q = Q, w.new = w.new, sw.new = sqrt(w.new), mu.new = mu.new,
             z.new = z.new, zw.new = z.new * sqrt(w.new), b.new = b.new,
             c.new = c.new, optlambda = optlambda, conv = TRUE)

  rm(K)
  rm(U)
  rm(Q)
  rm(Uv)
  rm(Qv)
  
  return(out)
}


cv.nng.subset = function(model, K, y, nbasis, basis.id, mscale, lambda0, lambda_theta, gamma, cv, nfold, one.std, obj)
{
  message("-- theta-step -- \n")
  message("proceeding... \n")
  n = length(y)
  d = length(mscale)
  
  Gw <- matrix(0, n, d)
  for (j in 1:d) {
    Gw[, j] = ((model$Uv[, , j] * sqrt(model$w.new)) %*% model$c.new) * (mscale[j]^(-2))
  }
  
  G <- matrix(0, n, d)
  for (j in 1:d) {
    G[, j] = (model$Uv[, , j] %*% model$c.new) * (mscale[j]^(-2))
  }
  
  uw = model$zw.new - model$sw.new
  
  h = rep(0, d)
  for (j in 1:d) {
    h[j] = n * lambda0 * ((t(model$c.new) %*% model$Uv[basis.id, , j]) %*% model$c.new)
  }
  
  # cv = "GCV" -> use only training data
  if(cv == "GCV"){
    init.theta = rep(1, d)
    len = length(lambda_theta)
    measure = numeric(NA, len)
    
    for (k in 1:len) {
      theta.new = .Call("wls_theta_step", Gw, uw, h/2, n, d, init.theta, n * lambda_theta[k] * gamma / 2, n * lambda_theta[k] * (1-gamma), PACKAGE = "cossonet")
      theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
      
      te_U = wsGram(Uv, theta.adj/mscale^2)
      testf = c(U %*% model$c.new + model$b.new)
      
      err = n * sum(model$w.new * (y - testf)^2)
      inv.mat = ginv(t(U) %*% U + lambda_theta[k] * model$Q)
      df = sum(diag(U %*% inv.mat %*% t(U)))
      measure[k] = err / (te_n - df)^2
    }
    
    # smoothing parameter selection
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
    min_id = which.min(measure_mean)
    optlambda = cand.lambda[min_id]
    
    if(show){
      plot(log(cand.lambda), measure, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
           ylim = range(measure), pch = 15, col = 'red', type = "b")
      abline(v = log(optlambda), lty = 2, col = "darkgray")
    }
  }
  
  
  if(cv == "mse"){
    
    # cross-validation
    init.theta = rep(1, d)
    len = length(lambda_theta)
    measure <- matrix(NA, nfold, len)
    fold = cvsplitID(n, nfold, y, family = obj$family)
    
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
        theta.new = .Call("wls_theta_step", Gw[tr_id,], uw[tr_id], h/2, tr_n, d, init.theta, tr_n * lambda_theta[k] * gamma / 2, tr_n * lambda_theta[k] * (1-gamma), PACKAGE = "cossonet")
        theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
        
        te_U = wsGram(te_Uv, theta.adj/mscale^2)
        testf = c(te_U %*% model$c.new + model$b.new)
        testmu = obj$linkinv(testf)
        
        if(obj$family == "gaussian") measure[fid, k] <- mean((testf - y[te_id])^2)
        if(obj$family == "binomial") measure[fid, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
        if(obj$family == "poisson") measure[fid, k] <- mean((y[te_id] - testf)^2)
        
        # if(cv == "mse"){
        #   testmu = obj$linkinv(testf)
        #
        #   if(obj$family == "gaussian") measure[fid, k] <- mean((testf - y[te_id])^2)
        #   if(obj$family == "binomial") measure[fid, k] <- mean(y[te_id] != ifelse(testmu < 0.5, 0, 1))
        #   if(obj$family == "poisson") measure[fid, k] <- mean((y[te_id] - testf)^2)
        # }
        
        # if(cv == "KL"){
        #   true_mu = obj$linkinv(f[te_id])
        #   measure[fid, k] <- KL(testf, true_mu, obj)
        # }
        
      }
    }
    
    rm(te_Uv)
    rm(te_U)
    
    
    # smoothing parameter selection
    if(obj$family == 'gaussian'){
      main = "Gaussian Family"
    }
    if(obj$family == 'binomial'){
      main = "Binomial Family"
    }
    if(obj$family == 'poisson'){
      main = "Poisson Family"
    }
    
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
    
    
    plot(log(lambda_theta), measure_mean, main = main, xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab,
         ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
    arrows(x0 = log(lambda_theta), y0 = measure_mean - measure_se,
           x1 = log(lambda_theta), y1 = measure_mean + measure_se,
           angle = 90, code = 3, length = 0.1, col = "darkgray")
    abline(v = log(optlambda), lty = 2, col = "darkgray")
    
  }
  
  theta.new = .Call("wls_theta_step", Gw, uw, h/2, n, d, init.theta, n * optlambda * gamma / 2, n * optlambda * (1-gamma), PACKAGE = "cossonet")
  theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
  
  f.new =  c(wsGram(model$Uv, theta.adj/mscale^2) %*% model$c.new + model$b.new)
  mu.new = obj$linkinv(f.new)
  
  if(obj$family == "gaussian") m = mean((f.new - y)^2)
  if(obj$family == "binomial") m <- mean(y != ifelse(mu.new < 0.5, 0, 1))
  if(obj$family == "poisson") m <- mean((y - f.new)^2)
  
  message("mse:", round(m, 4), "\n\n")
  
  out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.new)
  
  rm(G)
  rm(Gw)
  
  return(out)
}