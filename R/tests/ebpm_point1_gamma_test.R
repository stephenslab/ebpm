rm(list = ls())
source("../ebpm_point1_gamma.R")
source("../ebpm_gamma_mixture.R")
source("../point1_gamma.R")
source("../point_gamma.R")

source("../gammamix.R")

source("../utils.R")
library(mixsqp)

set.seed(123)
n = 10000
pi0 = 0.8
shape = 10
scale = 1

s = replicate(n, 2)
lam = replicate(n, 1)
n1 = round((1 - pi0)*n)
lam[1:n1] = rgamma(n = n1, shape = shape, scale = scale)
x = rpois(n = n, lambda = s * lam)



km = kmeans(x/s, 2)
c.id = which.max(km$centers)
pi0 = mean(km$cluster != c.id)
idx = which(km$cluster == c.id) ## indices of those non-one samples
c1 = mean(s[idx] * x[idx])
c2 = mean((s[idx] * x[idx])^2) - mean(s[idx]^2 * x[idx])
shape = c1^2/(c2 - c1^2)
scale = c1/(c2 - c1^2)

c(pi0, shape, scale)


start = proc.time()
fit = ebpm_point1_gamma(x = x, s = s,
                        fix_g = FALSE)
proc.time() - start


g_init = scale2gammamix_init(shape = c(seq(0.1, 0.9, 0.1), seq(1, 10, 1)), 
                             scale = c(seq(0.1, 0.9, 0.1), seq(1, 10, 1)))
start = proc.time()
fit = ebpm::ebpm_gamma_mixture(x = x, s = s,
                               g_init = g_init, 
                               fix_g = FALSE)
proc.time() - start

start = proc.time()
fit = ebpm::ebpm_point_gamma(x = x, s = s,
                         fix_g = FALSE)
proc.time() - start
# par(mfrow = c(2,1))
# plot(fit$posterior$mean)
# plot(fit$posterior$pi0_hat)
# fit$fitted_g
# fit$log_likelihood
