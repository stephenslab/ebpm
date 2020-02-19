rm(list = ls())
source("mle_two_gamma.R")
source("w_nb_mle.R")
library(stats)

obj.w.nb <- function(param, x, w){
  r = param[1]
  p = param[2]
  expr2 <- exp(r)
  loglik <- Rfast::Lgamma(x + expr2) - Rfast::Lgamma(x + 1) - lgamma(expr2) + expr2 * log(p) + x * log(1-p)
  loglik <- sum( w * loglik)
  return(-loglik)
}

nlm.w.nb <- function(x,w, r, p){
  param = c(r,p)
  start = proc.time()
  fit = nlm(f = obj.w.nb, p = c(r,p), w = w, x = x)
  runtime = proc.time() - start
  param <- c( fit$estimate[2], exp(fit$estimate[1]))
  names(param) <- c("success probability", "number of failures")
  list(iters = fit$iterations, loglik = - fit$minimum, param = param, runtime = runtime[[3]])
}






n = 10000
r = log(10); p = 0.1; w = replicate(n, 1/n)
x = rnbinom(n = n, size = exp(r), prob = p)

#fit0 = nlm.w.nb(x = x, w = w, r = r, p = p)
fit1 = w.negbin.mle(x = x, w = w, expr1 = exp(r),tol = 1e-3)
fit2 = w.negbin.mle2(x = x, w = w, expr1 = exp(r),tol = 1e-3)

#print(fit0)
print(fit1)
print(fit2)