## test
context("test_ebpm_point_gamma")
sim_spike_one <- function(pi, a, b){
  if(rbinom(1,1, pi)){return(0)}
  else{return(rgamma(1,shape = a, rate = b))}
}

simulate_pm <- function(s, param){
  pi = param[1]
  a = param[2]
  b  = param[3]
  lam = replicate(length(s), sim_spike_one(pi, a, b))
  x = rpois(length(s), s*lam)
  ll = -pg_nlm_fn(transform_param(param), x, s)
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
param =  c(pi, a, b)
sim = simulate_pm(s, param)

fit <- ebpm_point_gamma(sim$x, sim$s)

test_that("fitted loglikelihood > simulated  loglikelihood", {
  expect_gt(fit$log_likelihood, sim$log_likelihood)
})

test_that("RMSE: posterior  > MLE", {
  expect_lt(rmse(fit$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true))
  ## expect_gt(rmse(fit$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true)) ## this should  give error
})

## to do:
## add  comparison  with `ashr_pois` when it is updated.  




