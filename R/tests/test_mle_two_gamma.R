rm(list = ls())
library(ebpm)
source("two_gamma.R")
source("mle_two_gamma.R")
set.seed(123)
## test mle two gamma
simulate_tg_poisson <- function(s, pi0, shape1, scale1, shape2, scale2, 
								n = 1000, seed = 123){
	set.seed(seed)
	idx = rbinom(n = n, size = 1, prob = pi0)
	lam = idx * rgamma(n = n, shape = shape1, scale = scale1)
	+ (1 - idx) * rgamma(n = n, shape = shape2, scale = scale2)
	x = rpois(n = n, lambda = s*lam)
	return(list(x = x, lam = lam))
}

n = 10000
s = 1
pi0 = 0.7
shape1 = 0.5
scale1 = 1/10
shape2 = 5
scale2 = 2

maxiter = 30


g_init = list(pi0 = pi0, shape1 = shape1, scale1 = scale1, 
	shape2 = shape2, scale2 = scale2)

sim = simulate_tg_poisson(s, pi0, shape1, scale1, shape2, scale2, n)
x = sim$x
lam = sim$lam




start = proc.time() 
fit_tg = ebpm_two_gamma(x, s, g_init = g_init, n_iter = maxiter)
runtime_tg = proc.time() - start

start = proc.time() 
fit_tg_fast = mle_two_gamma(x, s, g_init, maxiter = maxiter, 
	tol_in = 1e-09, verbose = TRUE, get_progress = TRUE)
runtime_tg_fast = proc.time() - start

print("compare runtime")
print("tg")
print(runtime_tg)
print("tg_fast")
print(runtime_tg_fast)

print("compare loglikelihood")
print("tg")
print(fit_tg$log_likelihood)
print("tg_fast")
print(fit_tg_fast$progress[length(fit_tg_fast$progress)])


# [1] "compare runtime"
# [1] "tg"
#    user  system elapsed 
#  17.713   0.000  17.712 
# [1] "tg_fast"
#    user  system elapsed 
#   1.401   0.004   1.404 
# [1] "compare loglikelihood"
# [1] "tg"
# [1] -1498.943
# [1] "tg_fast"
# [1] -1498.943
