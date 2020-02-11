## test
rm(list = ls())
devtools::load_all(".")
context("test_ebpm_gh_gamma")
#library(ebpm)
set.seed(123)

library(gsl)
library(numDeriv)
## simulate data

## I make y not integers
n_z = 100
n_nz = 200
y = c(0.1 * runif(n_z), rpois(n = n_nz, lambda = 10) + runif(n_nz))
theta = c(replicate(n_z, 0.05), replicate(n_nz, 10.5))

#g_init = list(a = 0.5, b = 0.5, alpha = 0.01, phi = 0.001, gam = 0.01)
#g_init = NULL0
# fix_g = c(TRUE, TRUE, TRUE, TRUE, TRUE)
# fix_g = c(TRUE, TRUE, FALSE, FALSE, FALSE)

g_init = list(a = 0.1, b = 0.001, alpha = NULL, phi = NULL, gam = NULL)
fix_g = TRUE

fit = ebpm_gh_gamma(x = y, g_init = g_init, fix_g = fix_g)

fit2 = ebpm_two_gamma(x = y)
fit3 = ebpm_point_gamma(x = y)

#plot(y, fit$posterior$mean)
ix = 1:n_z

plot(y[1:n_z], fit2$posterior$mean[1:n_z])
abline(a = 0, b = 1, col = "blue")

plot(y[1:n_z], fit$posterior$mean[1:n_z])
abline(a = 0, b = 1, col = "red")



#plot(y[1:n_z], fit2$posterior$mean[1:n_z])

#plot(theta, fit2$posterior$mean)
# plot(y, fit3$posterior$mean)


fit$log_likelihood
fit2$log_likelihood
fit3$log_likelihood
fit$fitted_g
