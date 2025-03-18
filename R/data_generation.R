#' The function data_generation generates an example dataset for applying the cossonet function.
#'
#' @param n observation size.
#' @param p dimension.
#' @param rho a positive integer indicating the correlation strength for the first four informative variables.
#' @param SNR signal-to-noise ratio.
#' @param response the type of the response variable.
#' @param interaction the type of the effect, main or two-way interaction.
#'
#' @return a list of explanatory variables, response variables, and true functions.
#'
#' @examples
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
#' @export
#'
data_generation = function(n, p, rho, SNR,
                           response = c("continuous", "binary", "count", "survival"),
                           interaction = TRUE){

  if(response == "binary"){
    f1 = function(t) 3 * t
    f2 = function(t) pi * sin(pi * t) * 2
    f3 = function(t) 8 * t^3
    f4 = function(t) 4/ (exp(1) - 1) * exp(t)
  } else{
    f1 = function(t) t
    f2 = function(t) (2 * t - 1)^2
    f3 = function(t) sin(2 * pi * t) / (2 - sin(2 * pi * t))
    f4 = function(t) 0.1*sin(2 * pi * t) + 0.2*cos(2 * pi * t) + 0.3*sin(2 * pi * t)^2 + 0.4*cos(2 * pi * t)^3 + 0.5*sin(2 * pi * t)^3
    f5 = function(t) sin(2 * pi * t) / (2 - sin(2 * pi * t)) + .5
    f6 = function(t) 0.1*sin(2 * pi * t) + 0.2*cos(2 * pi * t) + 0.3*sin(2 * pi * t)^2 + 0.4*cos(2 * pi * t)^3 + 0.5*sin(2 * pi * t)^3 + .5
  }

  if(missing(response))
    type = "continuous"
  response = match.arg(response)

  if(missing(n)) n = 200
  if(missing(p)) p = 10
  if(missing(rho)) rho = 0.8
  if(missing(SNR)) SNR = 8

  if(p <= 4) stop("dimension size should be larger than 4.")

  t = 2
  pp = 4
  x = matrix(0, n, pp)
  W = matrix(runif(n * pp), n, pp)
  U = runif(n)
  for(j in 1:pp){
    x[, j] = (W[, j] + t * U)/(1 + t)
  }

  # Set the outer margins
  # par(oma = c(0, 0, 0, 0))
  # Set the inner margin
  # par(mar = c(4, 4, 3, 1))
  # par(mfrow = c(1,4))
  # curve(f1, 0, 1)
  # curve(f2, 0, 1)
  # curve(f3, 0, 1)
  # curve(f4, 0, 1)
  # par(mfrow = c(1,1))

  if(response == "continuous"){
    if(interaction){
      V_sig = var(1 * f1(x[,1])) + var(1 * f2(x[,2])) + var(2 * f3(x[,3])) + var(3 * f4(x[,4]))
      sd = sqrt(V_sig / SNR)
      f = 1 * f1(x[,1]) + 1 * f2(x[,2]) + 2 * f3(x[,3]) + 3 * f4(x[,4]) + rnorm(n, 0, sd)

      x_nois = matrix(runif(n * (p-pp), 0, 1), n, (p-pp))
      x = cbind(x, x_nois)
      out = list(x = x, f = f, y = f)
    }
    else{
      V_sig = var(f1(x[,1])) + var(f2(x[,2])) + var(f3(x[,3])) + var(f4(x[,4])) +
        var(f1(x[, 3] * x[, 4])) + var(f2(x[, 1] + x[, 3])) + var(f3(x[, 1] * x[, 2]))
      sd = sqrt(V_sig / SNR)
      f = f1(x[,1]) + f2(x[,2]) + 2 * f3(x[,3]) + 3 * f4(x[,4]) +
        5 * f1(x[, 3] * x[, 4]) + 3 * f2((x[, 1] * x[, 3])) + 2 * f3(x[, 1] * x[, 2])

      x_nois = matrix(runif(n * (p-pp), 0, 1), n, (p-pp))
      x = cbind(x, x_nois)
      out = list(x = x, f = f, y = f)
    }
  }

  if(response == "binary"){
    V_sig = var(f1(x[,1])) + var(f2(x[,2])) + var(f3(x[,3])) + var(f4(x[,4]))
    sd = sqrt(V_sig / SNR)
    f = f1(x[,1]) + f2(x[,2]) + f3(x[,3]) + f4(x[,4]) - 11 + rnorm(n, 0, sd)

    x_nois = matrix(runif(n * (p-4), 0, 1), n, (p-4))
    x = cbind(x, x_nois)
    prob = exp(f)/(exp(f) + 1)
    y = rbinom(n, 1, prob)

    out = list(x = x, f = f, y = y)
  }

  if(response == "count"){
    V_sig = var(1 * f1(x[,1])) + var(1 * f2(x[,2])) + var(2 * f5(x[,3])) + var(3 * f6(x[,4]))
    sd = sqrt(V_sig / SNR)
    f = 1 * f1(x[,1]) + 1 * f2(x[,2]) + 2 * f5(x[,3]) + 3 * f6(x[,4]) + rnorm(n, 0, sd)
    # print(sd)
    x_nois = matrix(runif(n * (p-pp), 0, 1), n, (p-pp))
    x = cbind(x, x_nois)

    f = f / 3
    mu = exp(f)
    y = rpois(n, mu)

    out = list(x = x, f = f, y = y)
  }

  if(response == 'survival'){

    V_sig = var(1 * f1(x[,1])) + var(1 * f2(x[,2])) + var(2 * f3(x[,3])) + var(3 * f4(x[,4]))
    sd = sqrt(V_sig / SNR)
    f = 1 * f1(x[,1]) + 1 * f2(x[,2]) + 2 * f3(x[,3]) + 3 * f4(x[,4]) + rnorm(n, 0, sd)

    x_nois = matrix(runif(n * (p - pp), 0, 1), n, (p - pp))
    x = cbind(x, x_nois)
    surTime = rexp(n, exp(f))
    cenTime = rexp(n, exp(-f) * runif(1, 4, 6))

    y = cbind(time = apply(cbind(surTime, cenTime), 1, min), status = 1 * (surTime < cenTime))

    out = list(x = x, f = f, y = y)
  }
  return(out)
}
