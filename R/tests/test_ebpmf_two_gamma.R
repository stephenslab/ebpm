rm(list = ls())
library("NNLM")
library("ebpm")
#devtools::load_all("~/Desktop/git/ebpm")
library("ebpmf.alpha")
#devtools::load_all("../ebpmf.alpha")
source("~/Desktop/git/ebpmf_demo/code/misc.R")
# source("../two_gamma.R")
# source("../utils.R")
# source("../mle_two_gamma5.R")
# source("../ebpm_two_gamma_fast5.R")

set.seed(123)
# n = 99
# p = 300
# k= 4

dataset = readRDS("./dataset_block.Rds")
X = dataset$X
L = dataset$L
init = dataset$init


n = nrow(X)
p = ncol(X)
k = ncol(L)


maxiter = 1000
verbose = TRUE
tol = -1

### two-gamma (this is very very slow, around 60 tims slower than others)
start = proc.time()
fit_ebpmf_tg = ebpmf.alpha::ebpmf(X, K = k,pm_func = ebpm::ebpm_two_gamma,
                                  pm_control = list(n_iter = 10),
                                  init = init, maxiter = maxiter, verbose = verbose, tol = tol)
runtime = proc.time()- start
fit_ebpmf_tg$runtime = runtime[[3]]

saveRDS(list(fit = fit_ebpmf_tg, sessionInfo = sessionInfo()), "test_ebpmf_two_gamma.Rds")