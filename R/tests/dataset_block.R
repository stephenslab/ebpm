rm(list = ls())
set.seed(123)
library("ebpm")
library("NNLM")

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

### initialiazation
init_lf = NNLM::nnmf(A = X, k = k, loss = "mkl", method = "scd", max.iter = 20, verbose = FALSE)
init = list(qg = initialize_qg_from_LF(L0 = init_lf$W, F0 = t(init_lf$H)))

out = list(X = X, lambda = lambda, L = L, F = F, init = init)
saveRDS(out, "./dataset_block.Rds")
