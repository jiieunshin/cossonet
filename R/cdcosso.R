#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param x Input matrix with $n$ observations and $p$ dimension.
#' @param y Response variable. The type can be continuous (default), binary class, non-negative count, survival.
#' @param family A character string representing one of the built-in families. The value depends on the type of response variable, `family='gaussian'` for continuous, `family='binomial'` forj binary class, `family='poisson'` for non-negative count , and `family='Cox'` for survival.
#' @param kernel Kernel function which is used to convert the input data for training and predicting. The four types is provided, `linear` (default), `gaussian`, `poly`, and `spline`.
#' @param nbasis The number of basis.
#' @param basis.id The index of basis.
#' @param effect Effect of the component to be analyzed, `effect = "main"` for main effect (default), and `effect = "interaction"` for interaction.
#' @param kparam Kernel parameter values that is used in gaussian kernel and polynomial kernel.
#' @param lambda0 A vector of tuning parameter to select optimal smoothing parameter.
#' @param lambda_theta A vector of tuning parameter to select optimal tuning parameter.
#' @param gamma The elastic net mixing parameter between 0 and 1. When `gamma = 1`, the penalty is to be LASSO. When `gamma = 0`, the penalty is to be ridge penalty. The default is `gamma = 0.95`.
#' @param scale Boolean for whether to scale the input data to range between 0 and 1.
#'
#' @return A list containing information about the fitted model.
#' @export

# x = tr_x
# y = tr_y
# obj = binomial()
# gamma = 1
# type = "spline"
# one.std = TRUE
# scale = T
# wt = rep(1, ncol(x))
# kparam = 1
# nfolds =5
# algo = "CD"
# lambda0 = exp(seq(log(2^{-6}), log(2^{2}), length.out = 20))
# lambda_theta = exp(seq(log(2^{-6}), log(2^{2}), length.out = 20))
cdcosso = function (x,
                    y,
                    family = c("gaussian", "binomial", "poisson", "Cox"),
                    nbasis, basis.id,
                    kernel = c("linear", "gaussian", "poly", "spline"),
                    effect = c("main", "interaction"),
                    kparam = 1,
                    lambda0 = exp(seq(log(2^{-6}), log(2^{2}), length.out = 20)),
                    lambda_theta = exp(seq(log(2^{-6}), log(2^{2}), length.out = 20)),
                    gamma = 0.95, scale = TRUE)
{
  n = nrow(x)
  colnames(x) = NULL
  rownames(x) = NULL
  # if(class(x)[1] != "data.frame")
    # stop("A input x must be matrix")

  # family
  family = match.arg(family)
  if(family == "gaussian")
    obj = gaussian()
  if(family == "binomial")
    obj =  binomial()
  if(family == "poisson")
    obj = poisson()

  if(missing(kernel))
    type = 'spline'
  else
    type = match.arg(kernel)

  if(missing(effect))
    effect = 'main'
  else
    effect = match.arg(kernel)

  if(effect == "interaction") kernel = paste0(kernel, "2")

  # if(missing(lambda0))
  #   lambda0 = exp(seq(log(2^{-11}), log(2^{2}), length.out = 20))
  #
  # if(missing(lambda_theta))
  #   lambda_theta = exp(seq(log(2^{-11}), log(2^{2}), length.out = 20))

  if (family == "Cox" & !all(match(c("time", "status"), dimnames(y)[[2]], 0))) {
    stop("Cox model requires a matrix with columns 'time' and 'status' as a response")
  }

  objnm = ifelse(family == 'gaussian' | family == 'binomial' | family == 'poisson', 'glm', family)

  wt = rep(1, ncol(x))

  # fitting
  out = switch(objnm,
               glm = cdcosso.glm(x, y, wt, nbasis, basis.id, lambda0, lambda_theta, gamma, obj, type, kparam, scale),
               Cox = cdcosso.cox(x, unlist(y[, "time"]), unlist(y[, "status"]), nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, kparam, scale)
               # Negbin, svm 추???
  )

  attr(out, "class") = "cdcosso"

  return(out)
}
