# ## test
# rm(list = ls())
# devtools::load_all(".")
# context("test_ebpm_two_gamma_fast3")
# library(ebpm)
# set.seed(123)

# simulate_two_gamma_poisson <- function(g_true, s, n_sample = 1000){
#   lam1 = rgamma(n = n_sample, shape = g_true$shape1, scale = g_true$scale1)
#   lam2 = rgamma(n = n_sample, shape = g_true$shape2, scale = g_true$scale2)
#   z = rbinom(n = n_sample, size = 1, prob = g_true$pi0)
#   lam = z * lam1 + (1-z) * lam2
#   x = rpois(n = n_sample, lambda = s * lam)
#   ll = compute_ll_two_gamma(x, s, pi1 = g_true$pi0, a1 = g_true$shape1, b1 = 1/g_true$scale1, 
#                             a2 = g_true$shape2, b2 = 1/g_true$scale2)
#   #return(list(x = x, lam = lam))
#   return(list(x = x, s= s, lam_true = lam, g_true = g_true, log_likelihood = ll))
# }

# rmse <- function(x,y){
#   return(sqrt(mean((x-y)^2)))
# }

# g_true = two_gamma(pi0 = 0.6, shape1 = 500, scale1 = 1/10, shape2 = 20, scale2 = 1/2)
# sim = simulate_two_gamma_poisson(g_true, s = 1, n_sample = 1000)
  
# fit_start_truth = ebpm::ebpm_two_gamma_fast3(sim$x, sim$s, g_init = g_true)

# fit <- ebpm::ebpm_two_gamma_fast3(sim$x, sim$s)


# fit_fixg <- ebpm::ebpm_two_gamma_fast3(sim$x, sim$s, g_init = fit$fitted_g, fix_g = T)


# test_that("fitted loglikelihood > simulated  loglikelihood", {
#   expect_gt(fit$log_likelihood, sim$log_likelihood)
# })

# test_that("RMSE: posterior  > MLE when fix pi0 at truth", {
#   expect_lt(rmse(fit_fixg$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true))
# })

# test_that("fitted loglikelihood > simulated  loglikelihood  when fix pi0 at truth", {
#   expect_gt(fit_fixg$log_likelihood, sim$log_likelihood)
# })

# test_that("RMSE: posterior  > MLE", {
#   expect_lt(rmse(fit$posterior$mean, sim$lam_true), rmse(sim$x/sim$s, sim$lam_true))
# })


# test_that("test fix_g", {
#   expect_false(any(fit$posterior != fit_fixg$posterior) || !all.equal(fit$log_likelihood,fit_fixg$log_likelihood))
# })


# print(fit$fitted_g)

# ## to do:
# ## add  comparison  with `ashr_pois` when it is updated.  





