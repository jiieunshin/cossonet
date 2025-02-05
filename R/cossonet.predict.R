#' The function `cossonet.predict` predicts predictive values for new data based on an object from the `cossonet` function.
#'
#' @param model The fitted cossonet object.
#' @param testx The new data set to be predicted.
#'
#' @return A list of predicted values for the new data set.
#' @export

cossonet.predict = function(model, testx)
{
  family = model$family
  if(family == "gaussian") obj = gaussian()
  if(family == "binomial") obj = binomial()
  if(family == "poisson") obj = poisson()

  nbasis = length(model$data$basis.id)
  te_n <- dim(testx)[1]

  if(class(testx)[1] == "data.frame") testx = matrix(unlist(testx), nrow = te_n)
  testx = apply(testx, 2, rescale)

  K = make_anovaKernel(testx, model$data$x[model$data$basis.id, ], model$data$kernel, model$data$kparam)
  d = K$numK

  R = array(NA, c(te_n, nbasis, d))
  for(j in 1:d){
    R[, , j] = K$K[[j]]
  }

  Rtheta <- combine_kernel(R, model$theta_step$theta.new/(model$data$wt^2))

  if(family != "Cox"){
    f.new = c(Rtheta %*% model$c_step$c.new + model$c_step$b.new)

    out = list(f.new = f.new, mu.new = obj$linkinv(f.new))
  }

  if(family == "Cox"){
    f.new = c(Rtheta %*% model$c_step$c.new)

    out = list(f.new = f.new)
  }

  rm(K)
  rm(R)
  rm(Rtheta)

  return(out)
}
