#' @title Empirical Bayes Poisson Mean with Mixture of Exponential as Prior
#' @describeIn estimate mixture weight and posterior mean
#' @import mixsqp
#' @export 
#'
## compute ebpm_exponential_mixture problem
ebpm_exponential_mixture <- function(x,s,m = 2, grid = NULL, seed = 123){
  set.seed(seed)
  if(is.null(grid)){grid <- select_grid_exponential(x,s,m)}
  b = grid$b
  a = grid$a
  tmp <-  compute_L(x,s,a, b)
  L =  tmp$L
  l_rowmax = tmp$l_rowmax
  fit <- mixsqp(L, control = list(verbose = F))
  ll_pi = sum(log(exp(l_rowmax) * L %*%  fit$x))
  pi = fit$x
  cpm = outer(x,a,  "+")/outer(s, b, "+")
  lam_pm = cpm %*% pi
  ll_lam = sum(dpois(x, s*lam_pm, log = T))
  return(list(pi = pi, lam_pm = lam_pm, ll_lam = ll_lam,ll_pi = ll_pi,L = L,grid = grid))
}


geom_seq <- function(low, up, m){
  N =  ceiling(log(up/low)/log(m)) + 1
  out  = low*m^(seq(1,N, by = 1)-1)
  return(out)
}

lin_seq <- function(low, up, m){
  out = seq(low, up, length.out = m)
  return(out)
}

## select grid for b_k
select_grid_exponential <- function(x, s, m = 2){
  ## mu_grid: mu =  1/b is the exponential mean
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min(xprime/s)
  mu_grid_max = 2*max(x/s)
  mu_grid = geom_seq(mu_grid_min, mu_grid_max, m)
  #mu_grid = lin_seq(mu_grid_min, mu_grid_max, m)
  b = 1/mu_grid
  a = rep(1, length(b))
  return(list(a= a, b = b))
}

## compute L matrix from data and selected grid
## L_ik = NB(x_i; a_k, b_k/b_k + s_i)
## but for computation in mixsqr, we can simplyfy it for numerical stability
compute_L <- function(x, s, a, b){
  prob = 1 - s/outer(s,b, "+")
  l = dnbinom(x,a,prob = prob, log = T) 
  l_rowmax  = apply(l,1,max)
  L = exp(l -  l_rowmax)
  return(list(L = L, l_rowmax = l_rowmax))
}