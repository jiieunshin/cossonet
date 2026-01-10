cv.getc.subset = function(K, time, status,  nbasis, basis.id, mscale, c.init, 
                          cand.lambda, type, cv, nfold, one.std, show)
{
  d = K$numK
  n <- length(status)
  len = length(cand.lambda)
  
  ## ---- Construct U, Q ----
  Uv <- array(NA, c(n, nbasis, d))
  for(j in 1:d) Uv[,,j] <- K$K[[j]][, basis.id]
  U <- combine_kernel(Uv, mscale)
  
  Qv <- array(NA, c(nbasis, nbasis, d))
  for(j in 1:d) Qv[,,j] <- K$K[[j]][basis.id, basis.id]
  Q <- combine_kernel(Qv, mscale)
  
  if(cv == "GCV"){
    
    ## ---- Stabilize Q ----
    EigQ <- eigen(Q)
    loop <- 0
    while(min(EigQ$values) < 1e-10 && loop < 10){
      Q <- Q + 1e-8 * diag(nbasis)
      EigQ <- eigen(Q)
      loop <- loop + 1
    }
    EigQ$values[EigQ$values < 0] = 1e-08
    pseudoX = U %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))
    
    measure = rep(0, len)
    c.init0 = c.init
    
    for (k in 1:len){
      c.init = c.init0
      
      if(is.null(c.init)){
        response <- survival::Surv(time = time, event = status)
        fit0 <- glmnet(pseudoX, response,
                       family = "cox",
                       alpha = 0,
                       standardize = TRUE)
        
        # c.init <- as.vector(coef(fit0, s = fit0$lambda[ceiling(length(fit0$lambda)*0.8)]))
        # c.init[!is.finite(c.init)] <- 0
        # c.init <- pmin(pmax(c.init, -1), 1)
        # lp0 <- as.vector(U %*% c.init)
        # cap_lp <- 3
        # maxlp <- max(abs(lp0), na.rm = TRUE)
        # if (is.finite(maxlp) && maxlp > cap_lp) {
        #   c.init <- c.init * (cap_lp / maxlp)
        # }
        c.init = as.vector(glmnet(pseudoX, response, family = "cox", lambda = 2^{-20}, alpha = 0, standardize = TRUE)$beta)
      }
      
      # eta = as.vector(exp(U %*% c.init))
      f.init <- as.vector(U %*% c.init)
      f.init <- pmin(pmax(f.init, -3), 3)         # clip threshold는 3~5 사이에서 실험 가능
      eta <- exp(f.init)
      coxgrad_results = coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
      w <- -attr(coxgrad_results, "diag_hessian")
      # z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
      
      # NaN/Inf/<=0 처리 + floor
      w[!is.finite(w) | w <= 0] <- 0
      w_floor <- 1e-6
      w[w < w_floor] <- w_floor
      g <- coxgrad_results
      g[!is.finite(g)] <- 0
      z <- f.init - g / w
      
      zw = z * sqrt(w)
      Uw = U * sqrt(w)
      sw = sqrt(w)
      
      fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, as.integer(n), as.integer(nbasis),
                  n * cand.lambda[k], PACKAGE = "cossonet")
      
      c.new = fit$c.new
      b.new = fit$b.new
      
      f.new = as.vector(U %*% c.new + b.new)
      f.new = pmin(pmax(f.new, -3), 3)
      eta.new <- exp(f.new)
      coxgrad.new = coxgrad(eta.new, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
      w.new <- -attr(coxgrad.new, "diag_hessian")
      w.new[!is.finite(w.new) | w.new <= 0] <- 0
      w.new[w.new < w_floor] <- w_floor
      
      Uw.new = U * sqrt(w.new)
      
      # err = n * sum(w.new * (time - f.new)^2)
      # inv.mat = ginv(t(U) %*% diag(w.new) %*% U + n * cand.lambda[k] * Q)
      # df = sum(diag(Uw.new %*% inv.mat %*% t(Uw.new)))
      # # measure[k] = err / (n - df)^2
      # denom <- (n - df)
      # if (!is.finite(denom) || denom <= 1e-8) denom <- 1e-8
      # measure[k] <- err / (denom^2)
      
      # err = n * sum(w.new * (time - f.new)^2)
      # inv.mat = ginv(t(U) %*% U + cand.lambda[k] * Q)
      # df = sum(diag(U %*% inv.mat %*% t(U)))
      # measure[k] = err / (n - df)^2
      
      # calculate ACV -----
      RS = RiskSet(time, status)
      GH = cosso:::gradient.Hessian.C(
        c.new, Uv, Qv,
        time, status, mscale, cand.lambda[k], RS)
      UHU <- U %*% cosso:::My_solve(GH$H, t(U))
      ACV_pen <- sum(status == 1) / n^2 *
        (sum(diag(UHU)) / (n - 1) - sum(UHU) / (n^2 - n))

      measure[k] = PartialLik(time, status, RS, U %*% c.new + b.new) + ACV_pen
      
    }
    
    optlambda <- cand.lambda[ which.min(measure) ]
    
    if(show){
      plot(log(cand.lambda), measure, type="b", pch=15, col = "red",
           xlab="log(lambda0)", ylab="GCV",
           main="")
      abline(v=log(optlambda), col="darkgray", lty=2)
    }
  }
  
  ## =========================================================
  ## 2) MSE Cross-validation option
  ## =========================================================
  if(cv == "mse" & nfold > 1){
  
    fold = cvsplitID(n, nfold, status, family = "binomial")
    measure = matrix(0, nfold, len)
    c.init0 <- c.init
    
    for(fid in 1:nfold){
      
      tr <- na.omit(as.vector(fold[, -fid]))
      te <- na.omit(as.vector(fold[, fid]))
      ntr <- length(tr); nte <- length(te)
      
      ## fold 안전성: event가 없으면 skip (PartialLik/ACV 불가)
      if (sum(status[tr] == 1) == 0 || sum(status[te] == 1) == 0) {
        measure[fid, ] <- NA_real_
        next
      }
      
      # train test split
      Utrv = array(NA, c(ntr, nbasis, d))
      for(j in 1:d){
        Utrv[, , j] = K$K[[j]][tr, basis.id]
      }
      
      Utr = combine_kernel(Utrv, mscale)
      
      Utev = array(NA, c(nte, nbasis, d))
      for(j in 1:d){
        Utev[, , j] = K$K[[j]][te, basis.id]
      }
      
      Ute <- combine_kernel(Utev, mscale)
      
      loop = 0
      EigQ = eigen(Q)
      while (min(eigen(Q)$values) < 0 & loop < 10) {
        loop = loop + 1
        Q = Q + 1e-08 * diag(nbasis)
        EigQ = eigen(Q)
      }
      if (loop == 10)
        EigQ$values[EigQ$values < 0] = 1e-08
      pseudoX = Utr %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))
      
      rm(Utrv)
      # rm(Ute)

      for (k in 1:len){
        response <- survival::Surv(time = time[tr], event = status[tr])
        c.init <- c.init0
        
        if(is.null(c.init)){
          fit0 <- glmnet(
            pseudoX, response,
            family = "cox",
            alpha = 0,
            standardize = TRUE
          )
          
          idx <- ceiling(length(fit0$lambda) * 0.8)
          lam_init <- fit0$lambda[idx]
          
          c.init <- as.vector(coef(fit0, s = lam_init))
          c.init[!is.finite(c.init)] <- 0
          c.init <- pmin(pmax(c.init, -1), 1)
          
          # fit0 <- glmnet(pseudoXtr, response,
          #                family = "cox",
          #                alpha = 0,
          #                standardize = TRUE)
          # 
          # c.init <- as.vector(coef(fit0, s = max(fit0$lambda)))
          # c.init <- pmin(pmax(c.init, -1), 1)
          # # c.init = as.vector(glmnet(pseudoXtr, response, family = "cox", lambda = cand.lambda[k], alpha = 1, standardize = FALSE)$beta)
        }
        
        # eta = as.vector(exp(Utr %*% c.init))
        f.init = as.vector(Utr %*% c.init)
        # f.init = pmin(pmax(f.init, -3), 3)
        eta = exp(f.init)
        
        coxgrad_results = coxgrad(eta, response, rep(1, ntr), 
                                  std.weights = FALSE, diag.hessian = TRUE
                                  )
        
        w = -attr(coxgrad_results, "diag_hessian")
        z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
        
        zw = z * sqrt(w)
        Uw = Utr * sqrt(w)
        sw = sqrt(w)
        
        fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, as.integer(ntr), as.integer(nbasis),
                    ntr * cand.lambda[k], PACKAGE = "cossonet")
        
        c.new = fit$c.new
        b.new = fit$b.new
        
        # calculate ACV -----  
        RSte = RiskSet(time[te], status[te])
        GHte = cosso:::gradient.Hessian.C(
          c.new, Utev, Qv,
          time[te], status[te],
          mscale, cand.lambda[k], RSte
        )
        UHU <- Ute %*% cosso:::My_solve(GHte$H, t(Ute))
        ACV_pen <- sum(status[te] == 1) / nte^2 *
          (sum(diag(UHU)) / (nte - 1) - sum(UHU) / (nte^2 - nte))
        
        measure[fid, k] <- PartialLik(time[te], status[te], RSte, Ute %*% c.new + b.new) + ACV_pen
        
        c.init <- NULL

        # f.te <- as.vector(Ute %*% c.use)
        # f.te = pmin(pmax(f.te, -3), 3)
        # RSte = cossonet:::RiskSet(time[te], status[te])
        # GHte = cosso::gradient.Hessian.C(c.new, Ute, Q, time[te], status[te], 
        #                                  mscale, cand.lambda[k], RSte)
        # UHU = Ute %*% cosso:::My_solve(GH$H, t(Ute))
        # ACV_pen = sum(status[te] == 1)/nte^2 * (sum(diag(UHU))/(nte - 1) - sum(UHU)/(nte^2 - nte))
        # measure[fid, k] =  PartialLik(time[te], status[te], RSte, Ute %*% c.new) + ACV_pen
      }
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
    
    
    mean_m <- colMeans(measure, na.rm=TRUE)
    se_m <- apply(measure,2,sd, na.rm=TRUE)/sqrt(nfold)
    
    id <- which.min(mean_m)
    if(one.std){
      cand <- which(mean_m <= mean_m[id] + se_m[id])
      cand <- cand[cand >= id]
      optlambda <- cand.lambda[max(cand)]
    } else {
      optlambda <- cand.lambda[id]
    }
    
    if(show){
      plot(log(cand.lambda), mean_m, pch=15, col = 'red',
           xlab="log(lambda)", ylab=paste0(nfold,"-CV"),
           ylim=range(c(mean_m-se_m,mean_m+se_m)))
      arrows(log(cand.lambda), mean_m-se_m,
             log(cand.lambda), mean_m+se_m,
             angle=90, code=3, length=0.08)
      abline(v=log(optlambda), col="darkgray", lty=2)
    }
    
    rm(Utr)
    rm(Ute)
    rm(RSte)
    rm(pseudoXtr)
  }
  
  ## ---- Final fit using optlambda ----
  if(nfold == 1) optlambda = cand.lambda
  
  response <- survival::Surv(time = time, event = status)
  
  if(is.null(c.init)){
    c.init = as.vector(glmnet(pseudoX, response, family = "cox", lambda = 2^{-20}, alpha = 0, standardize = TRUE)$beta)
  }
  
  RS = RiskSet(time, status)
  
  f.init <- as.vector(U %*% c.init)
  f.init <- pmin(pmax(f.init, -3), 3)         # clip threshold는 3~5 사이에서 실험 가능
  eta <- exp(f.init)
  coxgrad_results = coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
  w <- -attr(coxgrad_results, "diag_hessian")
  # z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
  
  # NaN/Inf/<=0 처리 + floor
  w[!is.finite(w) | w <= 0] <- 0
  w_floor <- 1e-6
  w[w < w_floor] <- w_floor
  g <- coxgrad_results
  g[!is.finite(g)] <- 0
  z <- f.init - g / w
  
  # eta = as.vector(exp(U %*% c.init))
  # coxgrad_results = coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
  # w = - attributes(coxgrad_results)$diag_hessian
  # z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
  
  zw = z * sqrt(w)
  Uw = U * sqrt(w)
  sw = sqrt(w)
  
  final = .Call("wls_c_step", zw, Uw, Q, c.init, sw, as.integer(n), as.integer(nbasis), 
                n * optlambda, PACKAGE = "cossonet")
  
  c.new = final$c.new
  b.new = final$b.new
  
  f.new = as.vector(U %*% c.new + b.new)
  coxgrad.new = coxgrad(exp(f.new), response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
  w.new = - attributes(coxgrad.new)$diag_hessian
  sw.new = sqrt(w.new)
  z.new = (exp(f.new) - 0) - ifelse(w.new != 0, -coxgrad.new/w.new, 0)
  zw.new = z.new * sqrt(w.new)
  
  GH.new = cosso::gradient.Hessian.C(c.new, Uv, Qv, time, status, 
                                     mscale, optlambda, RS)
  UHU = U %*% My_solve(GH.new$H, t(U))
  ACV_pen = sum(status == 1)/n^2 * (sum(diag(UHU))/(n - 1) - sum(UHU)/(n^2 - n))
  measure =  PartialLik(time, status, RS, U %*% c.new + b.new) + ACV_pen
  
  out = list(cv_error = measure, RS = RS, Uv = Uv, Q = Q, 
             w.new = w.new, sw.new = sw, 
             z.new = z.new, zw.new = zw.new,
             opt_lambda0 = optlambda, ACV_pen = ACV_pen,
             c.new = c.new, b.new = b.new
             )
  
  rm(Uv)
  rm(U)
  rm(Qv)
  rm(Q)
  rm(RS)
  rm(coxgrad_results)
  rm(GH.new)
  return(out)
}

cv.gettheta.subset = function (model, K, time, status, nbasis, basis.id, mscale,
                               lambda0, lambda_theta, gamma, cv, nfold, one.std)
{
  n = length(time)
  d = length(mscale)
  
  ## Precompute Uv for prediction
  Uv <- model$Uv
  Q <- model$Q
  
  ## Gw_j = w-scaled Gram direction
  Gw <- matrix(0, n, d)
  for(j in 1:d){
    Gw[,j] <- sqrt(model$w.new) * (Uv[,,j] %*% model$c.new) * (mscale[j]^(-2))
  }
  
  ## Working residual uw
  uw = model$zw.new - model$sw.new * model$b.new
  
  ## Penalty components h_j = c^T Q_j c
  h <- rep(0, d)
  for (j in 1:d) {
    Qj <- K$K[[j]][basis.id, basis.id]
    h[j] = as.vector(lambda0 * t(model$c.new) %*% Qj %*% model$c.new)
  }
  
  ## ---- CV ----
  if(cv == "GCV"){
    len = length(lambda_theta)
    measure = numeric(len)
    
    for (k in 1:len) {
      theta.new = .Call("wls_theta_step", Gw, uw, n * h/2, n, d, 
                        rep(1, d),
                        n * lambda_theta[k] * gamma / 2,
                        n * lambda_theta[k] * (1-gamma),
                        PACKAGE = "cossonet")
      
      U.new = wsGram(Uv, theta.new/mscale^2)
      f.new = as.vector(U.new %*% model$c.new + model$b.new)
      measure[k] =  cosso::PartialLik(time, status, RiskSet(time, status), f.new) + model$ACV_pen
      
      ## update U
      
      # rss_theta = sum((uw - Gw %*% theta.new )^2)
      # df = sum(diag(
      #   Gw %*%
      #     ginv(t(Gw)%*%Gw + diag(n * lambda_theta[k] * (1-gamma), d)) %*%
      #     t(Gw)
      # ))
      # measure[k] = n * rss_theta / (n - df)^2
      
      # testf = c(U %*% model$c.new + model$b.new
      # err = n * sum(model$w.new * (time - testf)^2)
      # inv.mat = ginv(t(U) %*% U + lambda_theta[k] * model$Q)
      # df = sum(diag(U %*% inv.mat %*% t(U)))
      # measure[k] = err / (n - df)^2
    }
    
    opt_lambda_theta <- lambda_theta[which.min(measure)]
    
    ylab <- expression("GCV(" * lambda[theta] * ")")
    plot(log(lambda_theta), measure,
         main = "Cox family",
         xlab = expression("log(" * lambda * ")"),
         ylab = ylab,
         type="b", pch=15, col="red")
    abline(v = log(opt_lambda_theta), col="darkgray", lty=2)
  }
  
  
  if(cv == "mse" & nfold == 1) stop("nfold should be >1.")
  
  if(cv == "mse" & nfold > 1){
    
    fold <- cvsplitID(n, nfold, family="gaussian")
    measure <- matrix(NA, nfold, len)
    
    for(fid in 1:nfold){
      tr <- na.omit(as.vector(fold[,-fid]))
      te <- na.omit(as.vector(fold[, fid]))
      ntr <- length(tr); nte <- length(te)
      
      Utev = array(NA, c(nte, nbasis, d))
      for(j in 1:d){
        Utev[, , j] = K$K[[j]][te, basis.id]
      }
      
      for (k in 1:len) {
        
        theta.new <- .Call(
          "wls_theta_step",
          Gw[tr,], uw[tr], ntr * h/2,
          ntr, d,
          rep(1,d),
          ntr*lambda_theta[k]*gamma/2,
          ntr*lambda_theta[k]*(1-gamma),
          PACKAGE="cossonet"
        )
        # response <- survival::Surv(time = time[tr_id], event = status[tr_id])
        # eta = exp(G[tr_id,] %*% init.theta)
        # coxgrad_results = coxgrad(eta, response, rep(1, tr_n), std.weights = FALSE, diag.hessian = TRUE)
        # w = - attributes(coxgrad_results)$diag_hessian
        # z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0) + lambda0 * G[tr_id,] %*% t(G[basis.id, ]) %*% model$c.new
        
        # theta.new = .Call("wls_theta_step", Gw[tr_id,], uw[tr_id], h/2, tr_n, d, init.theta, tr_n * lambda_theta[k] * gamma / 2, tr_n * lambda_theta[k] * (1-gamma), PACKAGE = "cossonet")
        # theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
        
        Ute.w = wsGram(Utev, theta.new/mscale^2)
        ftest = as.vector(Ute.w %*% model$c.new + model$b.new)
        measure[fid, k] =  cosso::PartialLik(time[te], status[te], RiskSet(time[te], status[te]), ftest) + model$ACV_pen
        
        # if(cv == "ACV") {
        #   ACV = cosso::PartialLik(time[tr_id], status[tr_id], RiskSet(time[tr_id], status[tr_id]), fhat) + model$ACV_pen
        #   measure[fid, k] = ACV
        # }
      }
    }
    
    ## ---- Select θ ----
    mean_m <- colMeans(measure, na.rm=TRUE)
    se_m <- apply(measure, 2, sd, na.rm=TRUE)/sqrt(nfold)
    min_id <- which.min(mean_m)
    
    if(one.std){
      cand_ids = which( mean_m <= mean_m[min_id] + se_m[min_id] )
      cand_ids = cand_ids[cand_ids >= min_id]
      std_id = max(cand_ids)
      opt_lambda_theta = lambda_theta[std_id]
    } else {
      opt_lambda_theta = lambda_theta[min_id]
    }
    
    plot(log(lambda_theta), mean_m, pch=15, ylim = range(mean_m-se_m, mean_m+se_m),
         xlab="log(theta)", ylab=paste0(nfold,"-CV"), col = "red")
    arrows(log(lambda_theta), mean_m-se_m,
           log(lambda_theta), mean_m+se_m,
           angle=90, code=3, length=0.08)
    abline(v=log(opt_lambda_theta), col="darkgray", lty=2)
    
  }
  
  # response = survival::Surv(time = time, event = status)
  # eta = exp(G %*% init.theta)
  # coxgrad_results <- coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
  # w = - attributes(coxgrad_results)$diag_hessian
  # z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0) + lambda0 * G %*% t(G[basis.id, ]) %*% model$c.new
  
  ## ---- Final theta ----
  theta.new <- .Call(
    "wls_theta_step",
    Gw, uw, h/2,
    n, d,
    rep(1,d),
    n*opt_lambda_theta*gamma/2,
    n*opt_lambda_theta*(1-gamma),
    PACKAGE="cossonet"
  )
  theta.adj <- ifelse(theta.new < 1e-6,0,theta.new)
  
  out = list(cv_error = measure, 
             optlambda_theta = opt_lambda_theta, 
             gamma = gamma, 
             theta.new = theta.adj
             )
  
  return(out)
}


# cv.getc.subset = function(K, time, status,  nbasis, basis.id, mscale, c.init, cand.lambda, type, cv, nfold, one.std, show)
# {
#   message("-- c-step -- \n")
#   message("proceeding... \n")
# 
#   d = K$numK
#   n <- length(status)
#   len = length(cand.lambda)
# 
#   Uv = array(NA, c(n, nbasis, d))
#   for(j in 1:d){
#     Uv[, , j] = K$K[[j]][, basis.id]
#   }
#   U <- combine_kernel(Uv, mscale)
# 
#   Qv = array(NA, c(nbasis, nbasis, d))
#   for(j in 1:d){
#     Qv[, , j] = K$K[[j]][basis.id, basis.id]
#   }
#   Q <- combine_kernel(Qv, mscale)
# 
#   # cv = "GCV" -> use only training data
#   if(cv == "GCV"){
#     pseudoX = U %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))
#     measure = numeric(len)
#     
#     for (k in 1:len){
#       response <- survival::Surv(time = time, event = status)
#       c.init = as.vector(glmnet(pseudoX, response, family = "cox", lambda = cand.lambda[k], alpha = 1, standardize = FALSE)$beta)
#       eta = exp(U %*% c.init)
#       coxgrad_results <- coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
#       w <- - attributes(coxgrad_results)$diag_hessian
#       z <- (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
#       
#       zw = z * sqrt(w)
#       Uw = U * w
#       sw = sqrt(w)
#       
#       fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, n, nbasis, n * cand.lambda[k], PACKAGE = "cossonet")
#       c.new = fit$c.new
#       
#       testf = c(U %*% c.new)
#       
#       err = n * sum(w * (time - testf)^2)
#       inv.mat = ginv(t(U) %*% U + cand.lambda[k] * Q)
#       df = sum(diag(U %*% inv.mat %*% t(U)))
#       measure[k] = err / (n - df)^2
#     }
#     
#     ylab = expression("GCV(" * lambda[0] * ")")
#     min_id = which.min(measure)
#     optlambda = cand.lambda[min_id]
#     
#     if(show){
#       plot(log(cand.lambda), measure, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
#            ylim = range(measure), pch = 15, col = 'red', type = "b")
#       abline(v = log(optlambda), lty = 2, col = "darkgray")
#     }
#     rm(inv.mat)
#   }
#   
#   
#   # cv = "mse" -> nfold
#   if(cv == "mse"){
#     
#     fold = cvsplitID(n, nfold, y, family = obj$family)
#     measure = matrix(NA, nfold, len)
#     for(fid in 1:nfold){
#       tr_id = as.vector(fold[, -fid])
#       te_id = fold[, fid]
#   
#       tr_id = tr_id[!is.na(tr_id)]
#       te_id = te_id[!is.na(te_id)]
#   
#       tr_n = length(tr_id)
#       te_n = length(te_id)
#   
#       tr_U = array(NA, c(tr_n, nbasis, d))
#       for(j in 1:d){
#         tr_U[, , j] = K$K[[j]][tr_id, basis.id]
#       }
#   
#       tr_U <- combine_kernel(tr_U, mscale)
#   
#       te_U = array(NA, c(te_n, nbasis, d))
#       for(j in 1:d){
#         te_U[, , j] = K$K[[j]][te_id, basis.id]
#       }
#   
#       te_U <- combine_kernel(te_U, mscale)
#   
#       # initialize
#       loop = 0
#       EigQ = eigen(Q)
#       while (min(eigen(Q)$values) < 0 & loop < 10) {
#         loop = loop + 1
#         Q = Q + 1e-08 * diag(nbasis)
#         EigQ = eigen(Q)
#       }
#       if (loop == 10)
#         EigQ$values[EigQ$values < 0] = 1e-08
#       pseudoX = tr_U %*% EigQ$vectors %*% diag(sqrt(1/EigQ$values))
#   
#       for (k in 1:len){
#         response <- survival::Surv(time = time[tr_id], event = status[tr_id])
#         c.init = as.vector(glmnet(pseudoX, response, family = "cox", lambda = cand.lambda[k], alpha = 1, standardize = FALSE)$beta)
#         eta = exp(tr_U %*% c.init)
#         coxgrad_results <- coxgrad(eta, response, rep(1, tr_n), std.weights = FALSE, diag.hessian = TRUE)
#         w <- - attributes(coxgrad_results)$diag_hessian
#         z <- (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
#   
#         zw = z * sqrt(w)
#         Uw = tr_U * w
#         sw = sqrt(w)
#   
#         fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, tr_n, nbasis, tr_n * cand.lambda[k], PACKAGE = "cossonet")
#         measure[fid, k] =  PartialLik(time[tr_id], status[tr_id], tr_RS, tr_U %*% fit$c.new)
#       }
#       # if(cv == "ACV"){
#       #   te_RS = RiskSet(time[te_id], status[te_id])
#       #   tr_RS = RiskSet(time[tr_id], status[tr_id])
#       #   test_GH = cosso::gradient.Hessian.C(fit$c.new, te_R, Q, time[te_id], status[te_id], mscale, cand.lambda[k], te_RS)
#       #   train_GH = cosso::gradient.Hessian.C(fit$c.new, tr_R, Q, time[tr_id], status[tr_id], mscale, cand.lambda[k], tr_RS)
#       #
#       #   UHU = tr_U %*% My_solve(train_GH$H, t(tr_U))
#       #   ACV_pen = sum(status[tr_id] == 1)/tr_n^2 * (sum(diag(UHU))/(tr_n - 1) - sum(UHU)/(tr_n^2 - tr_n))
#       #   measure[fid, k] = PartialLik(time[tr_id], status[tr_id], tr_RS, tr_U %*% fit$c.new) + ACV_pen
#       # }
#     }
#     
#     # smoothing parameter selection
#     measure_mean = colMeans(measure, na.rm = T)
#     measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(5)
#   
#     sel_id = which(!is.nan(measure_se) & measure_se != Inf)
#     measure_mean = measure_mean[sel_id]
#     measure_se = measure_se[sel_id]
#     cand.lambda = cand.lambda[sel_id]
#   
#     min_id = which.min(measure_mean)
#   
#     if(one.std){
#       cand_ids = which((measure_mean >= measure_mean[min_id]) &
#                          (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
#       cand_ids = cand_ids[cand_ids >= min_id]
#       std_id = max(cand_ids)
#       optlambda = cand.lambda[std_id]
#     } else{
#       optlambda = cand.lambda[min_id]
#     }
#   
#     ylab = expression("GCV(" * lambda[0] * ")")
#   
#     if(show){
#       plot(log(cand.lambda), measure_mean, main = "Cox family", xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
#            ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
#       arrows(x0 = log(cand.lambda), y0 = measure_mean - measure_se,
#              x1 = log(cand.lambda), y1 = measure_mean + measure_se,
#              angle = 90, code = 3, length = 0.1, col = "darkgray")
#       abline(v = log(optlambda), lty = 2, col = "darkgray")
#     }
#     
#     rm(tr_U)
#     rm(te_U)
#   }
#   
#   RS = RiskSet(time, status)
# 
#   response = survival::Surv(time = time, event = status)
#   c.init = as.vector(glmnet(U, response, family = "cox", lambda = optlambda, alpha = 1, standardize = FALSE)$beta)
#   eta = exp(U %*% c.init)
#   coxgrad_results <- coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
#   w = - attributes(coxgrad_results)$diag_hessian
#   z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0)
# 
#   zw = z * sqrt(w)
#   Uw = U * w
#   sw = sqrt(w)
# 
#   fit = .Call("wls_c_step", zw, Uw, Q, c.init, sw, n, nbasis, n * optlambda, PACKAGE = "cossonet")
# 
#   out = list(measure = measure, Uv = Uv, Q = Q, RS = RS, f.new = c(U %*% fit$c.new),
#              zw.new = zw, w.new = w, sw.new = sw, c.new = fit$c.new,
#              optlambda = optlambda)
# 
#   rm(K)
#   rm(Uv)
#   rm(U)
#   rm(Qv)
#   rm(Q)
#   rm(RS)
#   return(out)
# }
# 
# cv.gettheta.subset = function (model, K, time, status, nbasis, basis.id, mscale, lambda0, lambda_theta, gamma, cv, nfold, one.std)
#   {
#   message("-- theta-step -- \n")
#   message("proceeding... \n")
# 
#   n = length(time)
#   d = length(mscale)
# 
#   G <- matrix(0, n, d)
#   for (j in 1:d) {
#     G[, j] = (model$Uv[, , j] %*% model$c.new) * (mscale[j]^(-2))
#   }
# 
#   Gw <- matrix(0, n, d)
#   for (j in 1:d) {
#     Gw[, j] = ((model$Uv[, , j] * sqrt(model$w.new)) %*% model$c.new) * (mscale[j]^(-2))
#   }
# 
#   uw = model$zw.new - model$sw.new
# 
#   h = rep(0, d)
#   for (j in 1:d) {
#     h[j] = n * ((t(model$c.new) %*% model$Uv[basis.id, , j]) %*% model$c.new) * lambda0
#   }
#   
#   # cv = "GCV" -> use only training data
#   if(cv == "GCV"){
#     init.theta = rep(1, d)
#     len = length(lambda_theta)
#     measure = numeric(len)
#     
#     for (k in 1:len) {
#       theta.new = .Call("wls_theta_step", Gw, uw, h/2, n, d, init.theta, n * lambda_theta[k] * gamma / 2, n * lambda_theta[k] * (1-gamma), PACKAGE = "cossonet")
#       theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
#       
#       U = wsGram(model$Uv, theta.adj/mscale^2)
#       testf = c(U %*% model$c.new + model$b.new)
#       
#       err = n * sum(model$w.new * (y - testf)^2)
#       inv.mat = ginv(t(U) %*% U + lambda_theta[k] * model$Q)
#       df = sum(diag(U %*% inv.mat %*% t(U)))
#       measure[k] = err / (n - df)^2
#     }
#     
#     # smoothing parameter selection
#     if(obj$family == 'gaussian'){
#       main = "Gaussian Family"
#     }
#     if(obj$family == 'binomial'){
#       main = "Binomial Family"
#     }
#     if(obj$family == 'poisson'){
#       main = "Poisson Family"
#     }
#     
#     ylab = expression("GCV(" * lambda[0] * ")")
#     min_id = which.min(measure)
#     optlambda = lambda_theta[min_id]
#     
#     plot(log(lambda_theta), measure, main = main, xlab = expression("Log(" * lambda[0] * ")"), ylab = ylab,
#          ylim = range(measure), pch = 15, col = 'red', type = "b")
#     abline(v = log(optlambda), lty = 2, col = "darkgray")
#   }
#   
#   # cv = "mse" -> nfold
#   if(cv == "mse"){
#     
#     # cross-validation
#     init.theta = rep(1, d)
#     len = length(lambda_theta)
#     measure = matrix(NA, nfold, len)
#     fold = cvsplitID(n, nfold, time, family = "gaussian")
#   
#     for(fid in 1:nfold){
#       tr_id = as.vector(fold[, -fid])
#       te_id = fold[, fid]
#   
#       tr_id = tr_id[!is.na(tr_id)]
#       te_id = te_id[!is.na(te_id)]
#   
#       tr_n = length(tr_id)
#       te_n = length(te_id)
#   
#       te_Uv = array(NA, c(te_n, nbasis, d))
#       for(j in 1:d){
#         te_Uv[, , j] = K$K[[j]][te_id, basis.id]
#       }
#   
#       for (k in 1:len) {
#         response <- survival::Surv(time = time[tr_id], event = status[tr_id])
#         eta = exp(G[tr_id,] %*% init.theta)
#         coxgrad_results = coxgrad(eta, response, rep(1, tr_n), std.weights = FALSE, diag.hessian = TRUE)
#         w = - attributes(coxgrad_results)$diag_hessian
#         z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0) + lambda0 * G[tr_id,] %*% t(G[basis.id, ]) %*% model$c.new
#   
#         theta.new = .Call("wls_theta_step", Gw[tr_id,], uw[tr_id], h/2, tr_n, d, init.theta, tr_n * lambda_theta[k] * gamma / 2, tr_n * lambda_theta[k] * (1-gamma), PACKAGE = "cossonet")
#         theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
#   
#         te_U = wsGram(te_Uv, theta.adj/mscale^2)
#         fhat = c(te_U %*% model$c.new)
#         measure[fid, k] =  cosso::PartialLik(time[tr_id], status[tr_id], RiskSet(time[tr_id], status[tr_id]), fhat)
#   
#         # if(cv == "ACV") {
#         #   ACV = cosso::PartialLik(time[tr_id], status[tr_id], RiskSet(time[tr_id], status[tr_id]), fhat) + model$ACV_pen
#         #   measure[fid, k] = ACV
#         # }
#         }
#       }
#       
#       rm(te_Uv)
#       rm(te_U)
#       
#       # smoothing parameter selection
#       measure_mean = colMeans(measure, na.rm = T)
#       measure_se = apply(measure, 2, sd, na.rm = T) / sqrt(5)
#       
#       sel_id = which(!is.nan(measure_se) & measure_se != Inf)
#       measure_mean = measure_mean[sel_id]
#       measure_se = measure_se[sel_id]
#       lambda_theta = lambda_theta[sel_id]
#       
#       min_id = which.min(measure_mean)
#       
#       
#       if(one.std){
#         cand_ids = which((measure_mean >= measure_mean[min_id]) &
#                            (measure_mean <= (measure_mean[min_id] + measure_se[min_id])))
#         cand_ids = cand_ids[cand_ids >= min_id]
#         std_id = max(cand_ids)
#         optlambda = lambda_theta[std_id]
#       } else{
#         optlambda = lambda_theta[min_id]
#       }
#       
#       ylab = expression("GCV(" * lambda[theta] * ")")
#       
#       if(show){
#       plot(log(lambda_theta), measure_mean, main = "Cox family", xlab = expression("Log(" * lambda[theta] * ")"), ylab = ylab,
#            ylim = range(c(measure_mean - measure_se, measure_mean + measure_se)), pch = 15, col = 'red')
#       arrows(x0 = log(lambda_theta), y0 = measure_mean - measure_se,
#              x1 = log(lambda_theta), y1 = measure_mean + measure_se,
#              angle = 90, code = 3, length = 0.1, col = "darkgray")
#       abline(v = log(optlambda), lty = 2, col = "darkgray")
#       
#     }
#   }
#   
#   response = survival::Surv(time = time, event = status)
#   eta = exp(G %*% init.theta)
#   coxgrad_results <- coxgrad(eta, response, rep(1, n), std.weights = FALSE, diag.hessian = TRUE)
#   w = - attributes(coxgrad_results)$diag_hessian
#   z = (eta - 0) - ifelse(w != 0, -coxgrad_results/w, 0) + lambda0 * G %*% t(G[basis.id, ]) %*% model$c.new
# 
#   theta.new = .Call("wls_theta_step", Gw, uw, h/2, n, d, init.theta, n * optlambda * gamma / 2, n * optlambda * (1-gamma), PACKAGE = "cossonet")
#   theta.adj = ifelse(theta.new <= 1e-6, 0, theta.new)
# 
#   out = list(cv_error = measure, optlambda_theta = optlambda, gamma = gamma, theta.new = theta.adj)
# 
#   rm(G)
#   rm(Gw)
# 
#   return(out)
# }
# 
# soft_threshold = function(a, b){
#   return(ifelse(a > 0 & b < abs(a), a - b, 0))
# }
