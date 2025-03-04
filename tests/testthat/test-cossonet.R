library(testthat)
library(cossonet)

test_that("cossonet runs without errors", {
   data_generation(n = 50, p = 10, SNR = 9, response = "continuous")
})
