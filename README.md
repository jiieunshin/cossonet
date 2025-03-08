# cossonet

## Installation

We first load the library for `cossonet` and set a seed for reproducibility.

```{r, eval=FALSE, echo=FALSE, message=FALSE, warning=FALSE}
devtools::install_github("jiieunshin/cossonet")
library(cossonet)
set.seed(20250101)
```

## Data generation

The function `data_generation` generates example datasets with continuous response. We generate a training set with $n=200$ and $p=20$, and a test set with $n=1000$ and $p=20$.
```{r, eval=FALSE, echo=FALSE, message=FALSE, warning=FALSE}
tr = data_generation(n = 200, p = 20, SNR = 9, response = "continuous")
te = data_generation(n = 1000, p = 20, SNR = 9, response = "continuous")
```

## Model fitting

The function `cossonet` is the main function that fits the model. We have to input training set in this function. And Specific values are required to the arguments, such as `family`, `lambda0, and `lambda_theta`.

```{r, eval=FALSE, echo=FALSE, message=FALSE, warning=FALSE}
lambda0_seq = exp(seq(log(2^{-5}), log(2^{-1}), length.out = 20))
lambda_theta_seq = exp(seq(log(2^{-8}), log(2^{-5}), length.out = 20))

fit = cossonet(tr$x, tr$y, family = 'gaussian',
	       lambda0 = lambda0_seq,
	       lambda_theta = lambda_theta_seq
	       )
```

## Prediction

The function `cossonet.predict` is used to predict new data based on the fitted model. The output includes predicted values $\hat{f}$ (from `f.new`) and $\hat{\mu}$ (from `mu.new`) for the new data. The predicted value and predictive accuracy for the test set using our fitted model can be obtained by
```{r, eval=FALSE, echo=FALSE, message=FALSE, warning=FALSE}
pred = cossonet.predict(fit, te$x)
mean((te$f - pred$f.new)^2)
```
