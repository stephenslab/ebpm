## test
context("test_ebpm_two_gamma")
library(ebpm)
set.seed(123)
sim_two_gamma <- function(pi, a1, b1, a2, b2){
  if(rbinom(1,1, pi)){return(rgamma(1,shape = a1, rate = b1))}
  else{return(rgamma(1,shape = a2, rate = b2))}
}

simulate_pm <- function(s, param){
  param_trans = tg_transform_param(param)
  param = as.numeric(param)
  pi = param[1]
  a1 = param[2]
  b1  = 1/param[3]
  a2 = param[4]
  b2  = 1/param[5]
  lam = replicate(length(s), sim_two_gamma(pi, a1, b1, a2, b2))
  x = rpois(length(s), s*lam)
  ll = -tg_nlm_fn(param_trans, x, s)
  return(list(x = x, s= s, lam_true = lam, param = param, log_likelihood = ll))
}

rmse <- function(x,y){
  return(sqrt(mean((x-y)^2)))
}

n = 4000
s = replicate(n, 1)
pi  = 0.8

a1 = 0.01
b1  = 10

a2 = 100
b2 = 10
param =  list(pi0 = pi, shape1 = a1, scale1 = 1/b1, shape2 = a2, scale2 = 1/b2)

sim = simulate_pm(s, param)

fit <- ebpm::ebpm_two_gamma(sim$x, sim$s)

fit_fixpi <- ebpm::ebpm_two_gamma(sim$x, sim$s, pi0 = pi)


fit_fixg <- ebpm::ebpm_two_gamma(sim$x, sim$s, g_init = fit$fitted_g, fix_g = T)


test_that("fitted loglikelihood > simulated  loglikelihood", {
  expect_gt(fit$log_likelihood, sim$log_likelihood)
})

test_that("RMSE: posterior  > MLE when fix pi0 at truth", {
  expect_lt(rmse(fit_fixg$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true))
})

test_that("fitted loglikelihood > simulated  loglikelihood  when fix pi0 at truth", {
  expect_gt(fit_fixg$log_likelihood, sim$log_likelihood)
})

test_that("RMSE: posterior  > MLE", {
  expect_lt(rmse(fit$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true))
})


test_that("test fix_g", {
  expect_false(any(fit$posterior != fit_fixg$posterior) || fit$log_likelihood != fit_fixg$log_likelihood)
})




## to do:
## add  comparison  with `ashr_pois` when it is updated.  





