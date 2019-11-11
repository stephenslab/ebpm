## test
context("test_ebpm_point_gamma")
library(ebpm)
set.seed(123)
sim_spike_one <- function(pi, a, b){
  if(rbinom(1,1, pi)){return(0)}
  else{return(rgamma(1,shape = a, rate = b))}
}

simulate_pm <- function(s, param){
  param_trans = transform_param(param)
  param = as.numeric(param)
  pi = param[1]
  a = param[2]
  b  = 1/param[3]
  lam = replicate(length(s), sim_spike_one(pi, a, b))
  x = rpois(length(s), s*lam)
  ll = -pg_nlm_fn(param_trans, x, s)
  return(list(x = x, s= s, lam_true = lam, param = param, log_likelihood = ll))
}

rmse <- function(x,y){
  return(sqrt(mean((x-y)^2)))
}

n = 4000
s = replicate(n, 1)
pi  = 0.8
a = 100
b  = 1
param =  point_gamma(pi0 = pi, shape = a, scale = 1/b)

sim = simulate_pm(s, param)

fit <- ebpm::ebpm_point_gamma(sim$x, sim$s)

fit_fixpi <- ebpm::ebpm_point_gamma(sim$x, sim$s, pi0 = pi)


fit_fixg <- ebpm::ebpm_point_gamma(sim$x, sim$s, g_init = fit$fitted_g, fix_g = T)


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




