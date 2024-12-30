#' Load a matrix from a file
#'
#' The cdcosso function is a function that solves component selection using the Coordinate Descent algorithm.
#' This function can be applied to various response variables, continuous, count, binary, and survival analysis.
#' Because it is a type of nonparametric inference, various types of kernels can be selected.
#' To select hyperparameters, the function is designed to perform cross-validation.
#'
#' @param x Input matrix of $n$ by $p$, where each row is an observation. A matrix or dataframe. x must have at least two columns ($p>1$).
#' @param y The response variable. If family="gaussian" or family="poisson" (non-negative counts), it is quantitative. If family="binomial", it must be a vector with two levels. If family="cox", y must be a two-column matrix (or dataframe) with columns named 'time' and 'state'.
#' @param family The response type
#' @param wt The weights of the predictors. The default is `rep(1,ncol(x))`.
#' @param scale If TRUE, continuous predictors are rescaled to the interval `[0,1]`. The default is `TRUE`.
#' @param cv Measurement for cross-validation
#' @param nbasis The number of "knots" to choose from. If basis.id is provided, it is ignored.
#' @param basis.id An index that specifies the selected "knot".
#' @param kernel The kernel function. Four types are provided: `linear` (default), `gaussian`, `poly`, `spline`.
#' @param effect The effect of the component. `main` (default) for the main effect, `interaction` for interactions.
#' @param kparam Parameters for the kernel function. Used by Gaussian and polynomial kernels.
#' @param lambda0 A vector of lambda0 sequences. The default is `exp()`. This may need to be adjusted based on your data. Do not provide a single value for lambda0.
#' @param lambda_theta A vector of lambda_theta sequences. The default is `exp()`. This may need to be adjusted based on your data. Do not provide a single value for lambda_theta.
#' @param gamma Elastic mesh mixing parameter, `0 \leq \gamma \leq 1`. When `gamma = 1` it uses LASSO penalty, when `gamma = 0` it uses Ridge penalty. The default is `gamma = 0.95`.
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
                    wt = rep(1, ncol(x)),
                    scale = TRUE,
                    cv = NULL,
                    nbasis,
                    basis.id,
                    kernel = c("linear", "gaussian", "poly", "spline"),
                    effect = c("main", "interaction"),
                    kparam = 1,
                    lambda0 = exp(seq(log(2^{-6}), log(2^{2}), length.out = 20)),
                    lambda_theta = exp(seq(log(2^{-6}), log(2^{2}), length.out = 20)),
                    gamma = 0.95)
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

  if(is.null(cv) & family != "Cox")
    cv = "KL"

  if(is.null(cv) & family == "Cox")
    cv = "ACV"

  if(scale)
    x = apply(x, 2, rescale)

  if (family == "Cox" & !all(match(c("time", "status"), dimnames(y)[[2]], 0))) {
    stop("Cox model requires a matrix with columns 'time' and 'status' as a response")
  }

  objnm = ifelse(family == 'gaussian' | family == 'binomial' | family == 'poisson', 'glm', "Cox")

  # fitting
  out = switch(objnm,
               glm = cdcosso.glm(x, y, cv, wt, nbasis, basis.id, lambda0, lambda_theta, gamma, obj, type, kparam, scale),
               Cox = cdcosso.cox(x, unlist(y[, "time"]), unlist(y[, "status"]), cv, nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, kparam, scale)
  )

  attr(out, "class") = "cdcosso"

  return(out)
}
