#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param x Explanation variable matrix or data frame.
#' @param y Dependent variable vector or matrix or data frame containing time and status columns (for Cox model).
#' @param family Type of statistical model. Use one of the following strings: "gaussian", "binomial", "poisson", "negbin", "svm", or "Cox".
#' @param kernel Type of kernel function to use in case of SVM model. Use one of the following strings: "linear", "gaussian", "poly", "spline", "anova_gaussian", or "gaussian2".
#' @param algo Type of optimization algorithm. Use either the string "CD" (Coordinate Descent) or "QP".
#' @param wt Weight vector for each explanatory variable. The default is to use the same weight of 1 for all variables.
#' @param kparam Kernel parameter values to use for SVM models.
#' @param lambda0 Penalty parameter vector for Lasso and Ridge regression.
#' @param lambda_theta Vector of penalty parameters to be applied to different parts of the model.
#' @param M Vector of Lagrange multiplier.
#' @param gamma Gamma value used in Stochastic Search Optimization.
#' @param nfolds Number of folds for cross-validation.
#' @param one.std Logical value indicating whether to standardize explanatory variables.
#' @param scale A logical value indicating whether to scale the explanatory variable by min-max.
#' @param cpus Number of CPUs for parallel processing.
#'
#' @return A list containing information about the fitted model. Depending on the type of dependent variable, various information may be returned.
#' @export

# x = tr_x
# y = tr_y
# family = 'Cox'
# gamma = 0.8
# kernel = "gaussian"
# one.std = TRUE
# scale = T
# wt = rep(1, ncol(x))
# kparam = 1
# nfolds =5
cdcosso = function (x, y, family = c("gaussian", "binomial", "poisson", "negbin", "svm", "Cox"),
                    kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian", "gaussian2"),
                    algo = c("CD", "QP"), wt = rep(1, ncol(x)),
                    kparam = 1, lambda0, lambda_theta, M, gamma = 0.3, nfolds = 5, one.std = TRUE, scale = TRUE, cpus)
{
  n = nrow(x)
  colnames(x) = NULL
  rownames(x) = NULL
  if(class(x)[1] == "data.frame")
    x = as.matrix(x)

  # family
  family = match.arg(family)
  if(family == "gaussian")
    obj = gaussian()
  if(family == "binomial")
    obj =  binomial()
  if(family == "poisson")
    obj = poisson()
  if(family == "negbin"){
    link = poisson()$linkfun
    # if(missing(init.disp)){
    #   init.distp = NA
    # }
    obj = list(disp = NA, link = link)
  }


  if(missing(kernel))
    type = 'gaussian'
  else
    type = match.arg(kernel)

  if(missing(algo))
    algo = "CD"

  if(missing(lambda0)){
    lambda0 = exp(seq(log(2^{-40}), log(2^{10}), length.out = 40))
  }

  if(missing(lambda_theta))
    lambda_theta = exp(seq(log(2^{-40}), log(2^{4}), length.out = 40))

  if (scale){   # min-max scale
    x = apply(x, 2, rescale)
  }

  if (family == "Cox" & !all(match(c("time", "status"), dimnames(y)[[2]], 0))) {
    stop("Cox model requires a matrix with columns 'time' and 'status' as a response")
  }

  objnm = ifelse(family == 'gaussian' | family == 'binomial' | family == 'poisson', 'glm', family)

  # fitting
  out = switch(objnm,
               glm = cdcosso.glm(x, y, wt, lambda0, lambda_theta, gamma, obj, nfolds, one.std, type, kparam, algo)
               # Cox = cdcosso.cox(x, unlist(y[, "time"]), unlist(y[, "status"]), wt, lambda0, lambda_theta, gamma, nfolds, one.std, type, kparam, algo)
               # Negbin, svm 추가
  )
  return(out)
}
