rm(list = ls())
library("NNLM")
#library("ebpm")
devtools::load_all("~/Desktop/git/ebpm")
library("ebpmf.alpha")
#devtools::load_all("../ebpmf.alpha")
source("~/Desktop/git/ebpmf_demo/code/misc.R")

set.seed(123)
n = 99
p = 300
k= 4
mfac = 2.5 # controls PVE of dense factor
L = matrix(0, nrow=n, ncol=k)
F = matrix(0, nrow=p, ncol=k)
L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1
L[,4] = 1+mfac*runif(n)
F[1:(p/3),1] = 1+10*runif(p/3)
F[((p/3)+1):(2*p/3),2] = 1+10*runif(p/3)
F[((2*p/3)+1):p,3] = 1+10*runif(p/3)
F[,4]= 1+mfac*runif(p)
lambda = L %*% t(F)
X = matrix(rpois(n=length(lambda),lambda),nrow=n)

maxiter = 1000
verbose = TRUE
tol = -1

### initialiazation
init_lf = NNLM::nnmf(A = X, k = k, loss = "mkl", method = "scd", max.iter = 20, verbose = FALSE)
init = list(qg = initialize_qg_from_LF(L0 = init_lf$W, F0 = t(init_lf$H)))

### two-gamma (this is very very slow, around 60 tims slower than others)
start = proc.time()
fit_ebpmf_tg = ebpmf.alpha::ebpmf(X, K = k,pm_func = ebpm::ebpm_two_gamma_fast3,
                                  pm_control = list(n_iter = 10, control = list(tol_in = 1e-2)),
                                  init = init, maxiter = maxiter, verbose = verbose, tol = tol)
runtime = proc.time()- start
fit_ebpmf_tg$runtime = runtime[[3]]
