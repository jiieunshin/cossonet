#' Load a matrix from a file
#'
#' The cossonet function implements a nonparametric regression model that estimates nonlinear components.
#' This function can be applied to continuous, count, binary, and survival responses.
#' To use this function, the user must specify a family, kernel function, etc. For cross-validation, the sequence vectors `lambda0` and `lambda_theta` appropriate for the input data must also be specified.
#'
#' @param x Input matrix or data frame of $n$ by $p$. `x` must have at least two columns ($p>1$).
#' @param y A response vector with a continuous, binary, or count type. For survival responses, this should be a two-column matrix (or data frame) with columns called 'time' and 'status'.
#' @param family A distribution corresponding to the response type. `family="gaussian"` for continuous responses, `family="binomial"` for binary responses, `family="poisson"` for count responses, and `family="cox"` for survival responses.
#' @param wt The weights assigned to the explanatory variables. The default is `rep(1,ncol(x))`.
#' @param scale Boolean for whether to scale continuous explanatory variables to values between 0 and 1.
#' @param nbasis The number of "knots". If `basis.id` is provided, it is set to the length of `basis.id`.
#' @param basis.id The index of the "knot" to select.
#' @param kernel TThe kernel function. One of four types of `linear` (default), `gaussian`, `poly`, and `spline`.
#' @param effect The effect of the component. `main` (default) is the main effect, and `interaction` is the two-way interaction.
#' @param cv The criterion for selecting smoothing parameters.
#' @param nfold The number of folds to use in cross-validation is used to determine how many subsets to divide the data into for the training and validation sets.
#' @param kparam Parameters for Gaussian and polynomial kernel functions
#' @param lambda0 A vector of `lambda0` sequences. The default is a grid of 20 values `[2^{-10}, \dots, 2^{10}]` on an equally spaced logarithmic scale. This may need to be adjusted based on the input data. Do not set `\lambda0` as a single value.
#' @param lambda_theta A vector of `lambda` sequences. The default is a grid of 20 values `[2^{-10}, \dots, 2^{10}]` on an equally spaced logarithmic scale. This may need to be adjusted based on the input data. Do not set `lambda` as a single value.
#' @param gamma Elastic-net mixing parameter `0 \leq \gamma \leq 1`. If `gamma = 1`, the LASSO penalty is applied, and if `gamma = 0`, the Ridge penalty is applied. The default is `gamma = 0.95`.
#' @param one.std A logical value indicating whether to apply the "1-standard error rule." When set to `TRUE`, it applies to both the c-step and theta-step, selecting the simplest model within one standard error of the best model.

#' @return A list containing information about the fitted model.
#'
#' @examples
#' \donttest{
#' # Generate example data
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
#' }
#' @export
#'
cossonet = function (x,
                    y,
                    family = c("gaussian", "binomial", "poisson", "Cox"),
                    wt = rep(1, ncol(x)),
                    scale = TRUE,
                    nbasis,
                    basis.id,
                    kernel = c("linear", "gaussian", "poly", "spline"),
                    effect = c("main", "interaction"),
                    cv = c("GCV", "mse"),
                    nfold = 5,
                    kparam = 1,
                    lambda0 = exp(seq(log(2^{-29}), log(2^{-23}), length.out = 8)),
                    lambda_theta = exp(seq(log(2^{-16}), log(2^{-2}), length.out = 40)),
                    gamma = 0.95,
                    one.std = FALSE)
{
  n = nrow(x)
  colnames(x) = NULL
  rownames(x) = NULL
  if(!(class(x)[1] %in% c("matrix", "data.frame")))
    stop("A input x must be matrix or data frame.")

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

  if(missing(cv))
    cv = 'GCV'
  
  if(effect == "interaction") type = paste0(kernel, "2")

  if (family == "Cox" & !all(match(c("time", "status"), dimnames(y)[[2]], 0))) {
    stop("Cox model requires a matrix with columns 'time' and 'status' as a response")
  }

  objnm = ifelse(family == 'gaussian' | family == 'binomial' | family == 'poisson', 'glm', "Cox")

  # fitting
  out = switch(objnm,
               glm = cossonet.exp(x, y, wt, nbasis, basis.id, lambda0, lambda_theta, gamma, obj, type, cv, nfold, kparam, one.std, scale),
               Cox = cossonet.cox(x, unlist(y[, "time"]), unlist(y[, "status"]), nbasis, basis.id, wt, lambda0, lambda_theta, gamma, type, cv, nfold, kparam, one.std, scale)
  )

  attr(out, "class") = "cossonet"

  return(out)
}
