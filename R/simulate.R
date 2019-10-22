## functions  for simulate data for testing
#' @title Function for simulating from a mkixture of exponential
#' @import gtools

## simulate a poisson mean problem from mixture of exponential
simulate_pois_expmix  <-  function(n,seed = 123){
  set.seed(seed)
  ## simulate grid
  b = geom_seq(low = 0.1, up = 100, m = 2)
  #b = 1/geom_seq(low = 0.01, up = 10, m = 2)
  d = length(b)
  a = replicate(d,1)
  scale  = list(a = a, b = b)
  pi <- rdirichlet(1,rep(1/d, d))[1,]
  lam = replicate(n, sim_mgamma(a,b,pi))
  s = replicate(length(lam), 1)
  x  = rpois(length(lam),s*lam)
  tmp =  compute_L(x,s,a,b)
  L =  tmp$L
  l_rowmax = tmp$l_rowmax
  ll = sum(log(exp(l_rowmax) * L %*% matrix(pi, ncol = 1)))
  return(list(x =  x, s = s, lam = lam, g = gammamix(pi = pi, shape = scale$a, scale = 1/scale$b), ll = ll))
}


sim_mgamma <- function(a,b,pi){
  idx = which(rmultinom(1,1,pi) == 1)
  return(rgamma(1, shape = a[idx], rate =  b[idx]))
}

## simualte a poisson mean problem from point gamma

simulate_pois_point_gamma <- function(n, g_init = point_gamma(0.9,10,1/10), seed = 123){
  set.seed(seed)
  if(is.null(g_init)){
    g_init  = point_gamma(0.9,10,1/10)
  }
  lam = replicate(n, sim_pgamma(g_init$pi0, g_init$shape, g_init$scale))
  s = replicate(length(lam), 1)
  x  = rpois(length(lam),s*lam)
  ## compute ll
  ll =  ifelse(g_init$pi0 == 0, 
               -pg_nlm_fn_pi0(transform_param_pi0(g_init), x, s), 
               -pg_nlm_fn(transform_param(g_init), x, s))
  return(list(x =  x, s = s, lam = lam, g = g_init, ll = ll))
}

sim_pgamma <- function(pi0, shape, scale){
  if(rbinom(1,1,pi0) == 1){## means 0 for lambda
    return(0)
  }else{
    return(rgamma(n = 1, shape = shape, scale  = scale))
  }
}
  


