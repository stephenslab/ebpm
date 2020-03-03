rm(list = ls())
library(ebpm)
source("../two_gamma.R")
source("../mle_two_gamma4.R")
source("../utils.R")
set.seed(123)
## test mle two gamma
simulate_tg_poisson <- function(s, pi0, shape1, scale1, shape2, scale2, 
								n = 1000, seed = 123){
	set.seed(seed)
	idx = rbinom(n = n, size = 1, prob = pi0)
	lam = idx * rgamma(n = n, shape = shape1, scale = scale1) + 
	(1 - idx) * rgamma(n = n, shape = shape2, scale = scale2)
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

maxiter = 100


g_init = list(pi0 = 0.5, shape1 = shape1, scale1 = scale1, 
	shape2 = shape2, scale2 = scale2)

sim = simulate_tg_poisson(s, pi0, shape1, scale1, shape2, scale2, n)
x = sim$x
lam = sim$lam

control = list(nlm_setting = list(ndigit = 8, stepmax = 30, check.analyticals = FALSE),
						gradient = TRUE, hessian = FALSE)

print(control$gradient)
print(control$hessian)
start = proc.time() 
fit_tg_fast = mle_two_gamma(x, s, g_init, maxiter = maxiter, control = control,
	tol_in = 1e-09, verbose = FALSE, get_progress = TRUE)
runtime_tg_fast = proc.time() - start


# start = proc.time() 
# fit_tg = ebpm_two_gamma(x, s, g_init = g_init, n_iter = maxiter)
# runtime_tg = proc.time() - start



# print("compare runtime")
# print("tg")
# print(runtime_tg)
# print("tg_fast")
print(runtime_tg_fast)

# print("compare loglikelihood")
# print("tg")
# print(fit_tg$log_likelihood)
# print("tg_fast")
print(fit_tg_fast$progress[length(fit_tg_fast$progress)])


# > source("test_mle_two_gamma4.R")
# [1] FALSE
# [1] FALSE
#    user  system elapsed 
#   7.995   0.000   7.995 
# [1] -1498.654
# > source("test_mle_two_gamma4.R")
# [1] TRUE
# [1] FALSE
#    user  system elapsed 
#   6.049   0.000   6.050 
# [1] -1498.654
# > source("test_mle_two_gamma4.R")
# [1] TRUE
# [1] TRUE
#    user  system elapsed 
#   9.564   0.000   9.564 
# [1] -1498.654