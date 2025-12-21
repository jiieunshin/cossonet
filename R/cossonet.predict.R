#' The function `cossonet.predict` predicts predictive values for new data based on an object from the `cossonet` function.
#'
#' @param model The fitted cossonet object.
#' @param testx The new data set to be predicted.
#'
#' @return A list of predicted values for the new data set.
#'
#' @examples
#' \donttest{
#' set.seed(20250101)
#' tr = data_generation(n = 200, p = 20, SNR = 9, response = "continuous")
#' tr_x = tr$x
#' tr_y = tr$y
#'
#' te = data_generation(n = 1000, p = 20, SNR = 9, response = "continuous")
#' te_x = te$x
#' te_y = te$y
#'
#' # Fit the model
#' fit = cossonet(tr_x, tr_y, family = 'gaussian', gamma = 0.95, kernel = "spline", scale = TRUE,
#'       lambda0 = exp(seq(log(2^{-4}), log(2^{0}), length.out = 20)),
#'       lambda_theta = exp(seq(log(2^{-8}), log(2^{-6}), length.out = 20))
#'       )
#'
#' # Predict new dataset
#' pred = cossonet.predict(fit, te_x)
#' }
#' @export


cossonet.predict_v2 = function(model, testx){
  family = model$family
  obj <- switch(family,
                gaussian = gaussian(),
                binomial = binomial(),
                poisson = poisson(),
                Cox = NULL)
  
  # ensure matrix
  if(is.data.frame(testx)) testx <- as.matrix(testx)
  if(is.null(dim(testx)))  testx <- matrix(testx, nrow=1)
  
  te_n <- nrow(testx)
  
  # use training scaling!
  if(!is.null(model$data$Xmean) && !is.null(model$data$Xsd)){
    testx <- sweep(testx, 2, model$data$Xmean, "-")
    testx <- sweep(testx, 2, model$data$Xsd, "/")
  }
  
  nbasis = length(model$data$basis.id)
  
  K = make_anovaKernel(
    testx,
    model$data$x[model$data$basis.id, , drop=FALSE],
    model$data$kernel,
    model$data$kparam
  )
  
  d = K$numK
  R = array(NA, c(te_n, nbasis, d))
  for (j in 1:d) R[, , j] = K$K[[j]]
  
  Rtheta = combine_kernel(R, model$theta_step$theta.new/(model$data$wt^2))
  
  if(family!="Cox"){
    f.new = as.vector(Rtheta %*% model$c_step$c.new + model$c_step$b.new)
    return(list(f.new=f.new, mu.new=obj$linkinv(f.new)))
  }else{
    f.new = as.vector(Rtheta %*% model$c_step$c.new)
    return(list(f.new=f.new))
  }
}


# cossonet.predict = function(model, testx)
# {
#   family = model$family
#   if(family == "gaussian") obj = gaussian()
#   if(family == "binomial") obj = binomial()
#   if(family == "poisson") obj = poisson()
# 
#   nbasis = length(model$data$basis.id)
#   te_n <- dim(testx)[1]
# 
#   if(class(testx)[1] == "data.frame") testx = matrix(unlist(testx), nrow = te_n)
#   testx = apply(testx, 2, rescale)
# 
#   K = make_anovaKernel(testx, model$data$x[model$data$basis.id, ], model$data$kernel, model$data$kparam)
#   d = K$numK
# 
#   R = array(NA, c(te_n, nbasis, d))
#   for(j in 1:d){
#     R[, , j] = K$K[[j]]
#   }
# 
#   Rtheta <- combine_kernel(R, model$theta_step$theta.new/(model$data$wt^2))
# 
#   if(family != "Cox"){
#     f.new = c(Rtheta %*% model$c_step$c.new + model$c_step$b.new)
# 
#     out = list(f.new = f.new, mu.new = obj$linkinv(f.new))
#   }
# 
#   if(family == "Cox"){
#     f.new = c(Rtheta %*% model$c_step$c.new)
# 
#     out = list(f.new = f.new)
#   }
# 
#   rm(K)
#   rm(R)
#   rm(Rtheta)
# 
#   return(out)
# }
