## test
context("test_ebpm_gamma")
library(ebpm)
set.seed(123)

n = 1000
shape = 10
scale = 1/2
s = 1
lam = rgamma(n = n, shape = shape, scale = scale)
x = rpois(n = n, lambda = s * lam)

rmse <- function(lam, lamhat){
  sqrt(mean((lam - lamhat)^2))
}

fit <- ebpm::ebpm_gamma(x, s)

test_that("fitted loglikelihood > simulated  loglikelihood", {
  expect_gt(rmse(lam, x/s), rmse(lam, fit$posterior$mean))
})

