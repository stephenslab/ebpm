## test ebpm_two_gamma_fast3.R

rm(list = ls())
library(ebpm)
source("two_gamma.R")
source("utils.R")
source("mle_two_gamma.R")
source("ebpm_two_gamma_fast3.R")
set.seed(123)
## test mle two gamma
simulate_tg_poisson <- function(s, pi0, shape1, scale1, shape2, scale2,
								n = 1000, seed = 123){
	set.seed(seed)
	#browser()
	idx = rbinom(n = n, size = 1, prob = pi0)
	lam = idx * rgamma(n = n, shape = shape1, scale = scale1) + (1 - idx) * rgamma(n = n, shape = shape2, scale = scale2)
	x = rpois(n = n, lambda = s*lam)
	return(list(x = x, lam = lam))
}

n = 2000
s = 1
pi0 = 0.3
shape1 = 0.01
scale1 = 10
shape2 = 10
scale2 = 5

maxiter = 100

sim = simulate_tg_poisson(s, pi0, shape1, scale1, shape2, scale2, n)
x = sim$x
lam = sim$lam
# hist(lam)
hist(x)

start = proc.time()
fit_tg = ebpm_two_gamma(x, s, n_iter = maxiter, rel_tol = -1)
runtime_tg = proc.time() - start


control = list(tol_in = 1e-2)
start = proc.time()
fit_tg_fast = ebpm_two_gamma_fast3(x, s, n_iter = maxiter, control = control, verbose = FALSE)
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
print(fit_tg_fast$log_likelihood)

n_iter = min(length(fit_tg$progress), length(fit_tg_fast$progress))
print(n_iter)
plot(1:n_iter, fit_tg$progress[1:n_iter], col = "red")
lines(1:n_iter, fit_tg_fast$progress[1:n_iter], col = "blue")


# > source("test_ebpm_two_gamma_fast3.R")
# [1] "initialization time"
#    user  system elapsed 
#   0.084   0.000   0.085 
# [1] "mle fitting time"
#    user  system elapsed 
#   1.010   0.004   1.014 
# [1] "posterior computing time"
#    user  system elapsed 
#    0.03    0.00    0.03 
# [1] "compare runtime"
# [1] "tg"
#    user  system elapsed 
#  10.138   0.000  10.137 
# [1] "tg_fast"
#    user  system elapsed 
#   1.417   0.004   1.422 
# [1] "compare loglikelihood"
# [1] "tg"
# [1] -7246.649
# [1] "tg_fast"
# [1] -7246.649
# [1] 100
