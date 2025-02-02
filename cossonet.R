# Generated by `rjournal_pdf_article()` using `knitr::purl()`: do not edit by hand
# Please edit cossonet.Rmd to modify this file

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(plotly)
library(ggplot2)
library(palmerpenguins)
library(kableExtra)


## ----tab-comp-html, eval = knitr::is_html_output(), layout = "l-body-outset"----
# library(knitr)
# 
# table_data_html <- data.frame(
#   Response_Type = c("Gaussian", "Binomial", "Poisson", "Cox"),
#   cossonet = c("$\\checkmark$", "$\\checkmark$", "$\\checkmark$", "$\\checkmark$"),
#   cosso = c("$\\checkmark", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
#   LASSO = c("$\\checkmark", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
#   Ridge = c("$\\checkmark", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
#   Elastic_net = c("$\\checkmark$", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
#   main_effect = c("$\\checkmark$", "$\\checkmark$", "$\\checkmark$", "$\\checkmark$"),
#   two_way_interaction = c("$\\checkmark$", "$\\times$", "$\\checkmark$", "$\\times$")
# )
# 
# names(table_data_html) <- gsub("_", " ", names(table_data_html))
# 
# knitr::kable(table_data_html, format = "html", caption = "Comparison of the key capabilities of the cossonet and cosso packages.", align = "c")


## ----tab-comp-latex, eval = knitr::is_latex_output()--------------------------
library(knitr)
table_data_latex <- data.frame(
  Response_Type = c("Gaussian", "Binomial", "Poisson", "Cox"),
  cossonet = c("$\\checkmark$", "$\\checkmark$", "$\\checkmark$", "$\\checkmark$"),
  cosso = c("$\\checkmark", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
  LASSO = c("$\\checkmark", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
  Ridge = c("$\\checkmark", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
  Elastic_net = c("$\\checkmark$", "$\\checkmark$", "$\\times$", "$\\checkmark$"),
  main_effect = c("$\\checkmark$", "$\\checkmark$", "$\\checkmark$", "$\\checkmark$"),
  two_way_interaction = c("$\\checkmark$", "$\\times$", "$\\checkmark$", "$\\times$")
)

names(table_data_latex) <- gsub("_", " ", names(table_data_latex))

knitr::kable(table_data_latex, format = "latex", caption = "Comparison of the key capabilities of the cossonet and cosso packages.") %>% 
  kableExtra::kable_styling()


## ----call-library, echo = TRUE------------------------------------------------
devtools::install_github("jiieunshin/cossonet")
library(cossonet)
set.seed(20250101)


## ----gen-data, echo = TRUE----------------------------------------------------
tr = data_generation(n = 200, p = 20, SNR = 8, response = "regression",)
str(tr)

te = data_generation(n = 1000, p = 20, SNR = 8, response = "regression")
str(te)


## ----model-fitting, echo = TRUE, fig.weight = 15, fig.height = 3, fig.cap = "The five-fold CV plot for $\\lambda_0$ and $\\lambda$ from the `cossonet` run is shown. The horizontal axis represents the log of $\\lambda$, and the vertical axis represents the GCV. The red points are the average GCV for the validation set, while the vertical solid lines are the standard errors. If the 1-standard error rule is applied, the smoothing parameter is selected as the value corresponding to the dotted line in the current figure."----
fit = cossonet(tr$x, tr$y, family = "gaussian", 
              lambda0 = exp(seq(log(2^{-6}), log(2^{-3}), length.out = 20)),
              lambda_theta = exp(seq(log(2^{-8}), log(2^{-4}), length.out = 20))
              )


## ----tab-function-html, eval = knitr::is_html_output(), layout = "l-body-outset"----
# library(knitr)
# library(kableExtra)
# kable(
#   data.frame(
#     Argument = c("x", "y", "family", "wt", "scale", "cv", "nbasis", "basis.id", "kernel", "effect", "kparam", "lambda0", "lambda", "gamma"),
#     Description = c(
#       "Input matrix of size $n$ by $p$, where each row represents an observation. It can be a matrix or data frame. `x` must have at least two columns $(p \\gt 1)$.",
#       "The response variable. If `family=\"gaussian\"` or `family=\"poisson\"` (non-negative counts), it is quantitative. If `family=\"binomial\"`, `y` must be a vector with two levels. If `family=\"cox\"`, `y` must be a two-column matrix (or data frame) with columns named `time` and `state`.",
#       "The type of the response variable.",
#       "The weights of the predictors. The default is `rep(1, ncol(x))`.",
#       "If `TRUE`, continuous predictors are rescaled to the interval $[0, 1]$. The default is `TRUE`.",
#       "A measurement for cross-validation.",
#       "The number of \"knots\" to choose from. If `basis.id` is provided, it is ignored.",
#       "An index that specifies the selected \"knot\".",
#       "The kernel function. Four types are provided: `linear`, `gaussian`, `poly`, and `spline` (default). \n\n - Linear kernel: $K(x, y) = \\langle x, y\\rangle$ \n\n - Gaussian kernel: $K(x, y) = \\exp(-\\kappa\\langle x, y\\rangle)$ \n\n - Polynomial kernel: $K(x, y) = (\\langle x, y\\rangle +1)^\\kappa$ \n\n - Spline kernel: $K(x, y) = U(x, y)$",
#       "The effect of the component. `main` (default) for the main effect, `interaction` for two-way interactions.",
#       "Parameter $\\kappa$ for the kernel function. Used by Gaussian and polynomial kernels.",
#       "A vector of $\\lambda_0$ sequences. The default is a grid of 20 values $[2^{-10}, \\dots, 2^{10}]$ on an equally spaced logarithmic scale. This may need to be adjusted based on the input data. Do not provide a single value for $\\lambda_0$.",
#       "A vector of $\\lambda$ sequences. The default is a grid of 20 values $[2^{-10}, \\dots, 2^{10}]$ on an equally spaced logarithmic scale. This may need to be adjusted based on the input data. Do not provide a single value for $\\lambda$.",
#       "Elastic mesh mixing parameter, $0 \\leq \\gamma \\leq 1$. When `gamma=1`, it uses LASSO penalty, and when `gamma=0`, it uses ridge penalty. The default is `gamma=0.95`."
#     )
#   ),
#   caption = "Summary of the arguments of functions in `cossonet`.",
#   escape = FALSE
# ) %>%
#   kable_styling()


## ----tab-function-latex, eval = knitr::is_latex_output()----------------------
library(knitr)
kable(
  data.frame(
    Argument = c("x", "y", "family", "wt", "scale", "cv", "nbasis", "basis.id", "kernel", "effect", "kparam", "lambda0", "lambda", "gamma"),
    Description = c(
      "Input matrix of size $n$ by $p$, where each row represents an observation. It can be a matrix or data frame. \\texttt{x} must have at least two columns ($p>1$).",
      "The response variable. If \\texttt{family=\"gaussian\"} or \\texttt{family=\"poisson\"} (non-negative counts), it is quantitative. If \\texttt{family=\"binomial\"}, it must be a vector with two levels. If \\texttt{family=\"cox\"}, \\texttt{y} must be a two-column matrix (or data frame) with columns named 'time' and 'state'.",
      "The type of the response variable.",
      "The weights of the predictors. The default is \\texttt{rep(1, ncol(x))}.",
      "If \\texttt{TRUE}, continuous predictors are rescaled to the interval $[0, 1]$. The default is \\texttt{TRUE}.",
      "A measurement for cross-validation.",
      "The number of \"knots\" to choose from. If \\texttt{basis.id} is provided, it is ignored.",
      "An index that specifies the selected \"knot\".",
      "The kernel function. Four types are provided: \\texttt{linear}, \\texttt{gaussian}, \\texttt{poly}, and \\texttt{spline} (default). - Linear kernel: $K(x, y) = \\langle x, y \\rangle$ - Gaussian kernel: $K(x, y) = \\exp(-\\kappa \\langle x, y \\rangle)$ - Polynomial kernel: $K(x, y) = (\\langle x, y \\rangle +1)^{\\kappa}$ - Spline kernel: $K(x, y) = U(x, y)$",
      "The effect of the component. \\texttt{main} (default) for the main effect, \\texttt{interaction} for two-way interactions.",
      "Parameter $\\kappa$ for the kernel function. Used by Gaussian and polynomial kernels.",
      "A vector of $\\lambda_0$ sequences. The default is a grid of 20 values $\\mathtt{[2^{-10}, \\dots, 2^{10}]}$ on an equally spaced logarithmic scale. This may need to be adjusted based on the input data. Do not provide a single value for $\\lambda_0$.",
      "A vector of $\\lambda$ sequences. The default is a grid of 20 values $\\mathtt{[2^{-10}, \\dots, 2^{10}]}$ on an equally spaced logarithmic scale. This may need to be adjusted based on the input data. Do not provide a single value for $\\lambda$.",
      "Elastic mesh mixing parameter, $0 \\leq \\gamma \\leq 1$. When \\texttt{gamma=1}, it uses LASSO penalty, and when \\texttt{gamma=0}, it uses ridge penalty. The default is \\texttt{gamma=0.95}."
    )
  ),
  caption = "Summary of the arguments of functions in `cossonet`."
) %>%
  kable_styling()


## ----fit-str, echo = TRUE-----------------------------------------------------
str(fit)


## ----theta-new, echo = TRUE---------------------------------------------------
fit$theta_step$theta.new


## ----cossonet-predict, echo = TRUE--------------------------------------------
pred = cossonet.predict(fit, te$x)
str(pred)


## ----mse, echo = TRUE---------------------------------------------------------
mean((te$y - pred$f.new)^2)


## ----fig-nonlinear, echo = FALSE, fig.cap = "Nonlinear functions for generating simulation data."----
knitr::include_graphics("figures/nonlinear_funcs.png")


## ----simulate-table, echo=FALSE, message=FALSE, warning=FALSE-----------------
library(kableExtra)

data <- data.frame(
  n = rep(c(100, 200, 400), each = 8),
  p = rep(c(20, 40, 80, 160), 6),
  gamma = rep(c(0.95, 1), each = 4, times = 3),
  TP = c(3.25, 3.08, 2.57, 2.22, 2.69, 2.17, 3.16, 2.27, 
         3.69, 3.55, 2.68, 2.47, 2.71, 2.31, 3.15, 2.19,
         3.86, 3.86, 3.08, 2.98, 2.74, 2.33, 3.05, 2.34),
  FP = c(0.40, 0.35, 0.22, 0.22, 0.24, 0.22, 0.29, 0.21, 
         0.41, 0.34, 0.17, 0.17, 0.12, 0.07, 0.09, 0.03,
         0.19, 0.12, 0.05, 0.05, 0.02, 0.00, 0.00, 0.00),
  F1_score = c(0.8449, 0.8207, 0.7442, 0.6737, 0.7657, 0.6612, 0.8656, 0.6860, 
               0.9121, 0.8984, 0.7705, 0.7302, 0.7841, 0.7128, 0.8636, 0.6965,
               0.9594, 0.9669, 0.8546, 0.8390, 0.8008, 0.7258, 0.8548, 0.7297),
  MSE = c(2.6977, 2.6509, 3.0784, 3.0580, 2.9984, 2.9549, 3.0741, 3.0085,
          2.4456, 2.4333, 2.7630, 2.7469, 2.7027, 2.7296, 2.7004, 2.7403,
          2.3610, 2.3623, 2.5609, 2.5517, 2.6395, 2.6493, 2.6504, 2.7578),
  SE_TP = rep(c(0.0757, 0.0800, 0.0795, 0.0786, 0.0748, 0.0817, 0.0748, 0.0763), 3),
  SE_FP = rep(c(0.0865, 0.0796, 0.0613, 0.0561, 0.0553, 0.0596, 0.0656, 0.0556), 3),
  SE_F1 = rep(c(0.0143, 0.0139, 0.0154, 0.0166, 0.0134, 0.0173, 0.0136, 0.0158), 3),
  SE_MSE = rep(c(0.0293, 0.0282, 0.0389, 0.0376, 0.0364, 0.0332, 0.0413, 0.0364), 3)
)

data %>%
  kbl(
    caption = "Simulated results of continuous response for TP, FP, F1-score, and MSE for `cossonet` at γ = 1 and 0.95, with standard errors in parentheses.",
    digits = 4
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "center"
  ) %>%
  column_spec(3:10, bold = TRUE)

