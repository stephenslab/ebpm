## functions  for simulate data for testing
#' @title Function for simulating from a mkixture of exponential
#' @import gtools

## simulate a poisson mean problem
simulate_pois_expmix  <-  function(n,seed = 123){
  set.seed(seed)
  ## simulate grid
  b = geom_seq(low = 0.1, up = 100, m = 2)
  d = length(b)
  a = replicate(d,1)
  scale  = list(a = a, b = b)
  pi <- rdirichlet(1,rep(1/d, d))[1,]
  lam = replicate(n, sim_mgamma(a,b,pi))
  s = replicate(length(lam), 1)
  x  = rpois(length(lam),s*lam)
  ll_lam = sum(dpois(x, s*lam, log = T))
  tmp =  compute_L(x,s,a,b)
  L =  tmp$L
  l_rowmax = tmp$l_rowmax
  ll = sum(log(exp(l_rowmax) * L %*% matrix(pi, ncol = 1)))
  return(list(x =  x, s = s, lam = lam, g = gammamix(pi = pi, a = a, b = b), ll = ll))
}


sim_mgamma <- function(a,b,pi){
  idx = which(rmultinom(1,1,pi) == 1)
  return(rgamma(1, shape = a[idx], rate =  b[idx]))
}



