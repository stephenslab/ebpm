context("test_ebpm_exponential_mixture")
library(gtools)
## testing setup
sim_mgamma <- function(a,b,pi){
  idx = which(rmultinom(1,1,pi) == 1)
  return(rgamma(1, shape = a[idx], rate =  b[idx]))
}

## simulate a poisson mean problem
simulate_pm  <-  function(n, d, seed = 123){
  set.seed(seed)
  ## simulate grid
  a = replicate(d,1)
  b = 10*runif(d)
  grid  = list(a = a, b = b)
  pi <- rdirichlet(1,rep(1/d, d))
  lam_true = replicate(n, sim_mgamma(a,b,pi))
  s = replicate(length(lam_true), 1)
  x  = rpois(length(lam_true),s*lam_true)
  #ll_lam = sum(dpois(x, s*lam_true, log = T))
  #browser()
  tmp =  compute_L(x,s,a,b)
  L =  tmp$L
  l_rowmax = tmp$l_rowmax
  ll_pi = sum(log(exp(l_rowmax) * L %*% matrix(pi, ncol = 1)))
  return(list(x =  x, s = s, lam_true = lam_true, pi = pi, log_likelihood = ll_pi, grid = grid))
}

compute_L <- function(x, s, a, b){
  prob = 1 - s/outer(s,b, "+")
  l = dnbinom(x,a,prob = prob, log = T) 
  l_rowmax  = apply(l,1,max)
  L = exp(l -  l_rowmax)
  return(list(L = L, l_rowmax = l_rowmax))
}

rmse <- function(x,y){
  return(sqrt(mean((x-y)^2)))
}

##  simiulate data
n  = 400
d =  20
sim = simulate_pm(n, d)


## fit with ebpm_exponential_mixture
m = 1.1

fit       =  ebpm::ebpm_exponential_mixture(x = sim$x, s = sim$s, scale = "estimate", g_init = NULL, fix_g = F, m = 2, control = NULL)

fit_scale =  ebpm::ebpm_exponential_mixture(x = sim$x, s = sim$s, scale = fit$fitted_g$scale, g_init = NULL, fix_g = F, m = 2, control = NULL)

fit_init =  ebpm::ebpm_exponential_mixture(x = sim$x, s = sim$s, g_init = fit$fitted_g, fix_g = F, m = 2, control = NULL)

fit_fix =  ebpm::ebpm_exponential_mixture(x = sim$x, s = sim$s,  g_init = fit$fitted_g, fix_g = T, m = 2, control = NULL)


test_that("fitted loglikelihood > simulated  loglikelihood", {
  expect_gt(fit$log_likelihood, sim$log_likelihood)
})

test_that("RMSE: posterior  > MLE", {
  expect_lt(rmse(fit$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true))
  ## expect_gt(rmse(fit$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true)) ## this should  give error
})


test_that("test scale",{
  expect_false(any(fit$fitted_g$a != fit_scale$fitted_g$a)||
                 any(fit$fitted_g$b != fit_scale$fitted_g$b))
})

test_that("test init",{
  expect_false(any(fit$fitted_g$a != fit_init$fitted_g$a)||
                 any(fit$fitted_g$b != fit_init$fitted_g$b)||
                 fit$log_likelihood > fit_init$log_likelihood)
})

test_that("test fix g", {
  expect_false(any(fit$posterior != fit_fix$posterior) || 
                 !all.equal(fit$log_likelihood, fit_fix$log_likelihood,tolerance = 1e-5))
})


## to do:
## add  comparison  with `ashr_pois` when it is updated.  









