## description
## this solves ebpm problem, with mixture of exponential distribution as prior:
## g(.) = \sum_k \pi_k exp(.;b_k), where b_k is rate of exponential



library(mixsqp) ## for solving mixture problem
library(gtools) ## for using rdirichlet
library(ggplot2)

## generate a geometric sequence: x_n = low*m^{n-1} up to x_n < up
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
select_grid <- function(x, s, m = 2){
  ## mu_grid: mu =  1/b is the exponential mean
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min(xprime/s)
  mu_grid_max = 2*max(x/s)
  mu_grid = geom_seq(mu_grid_min, mu_grid_max, m)
  #mu_grid = lin_seq(mu_grid_min, mu_grid_max, m)
  b = 1/mu_grid
  return(list(b = b))
}

## compute L matrix from data and selected grid
## L_ik = (s_i^x_i* b_k)/(s_i + b_k)^(x_i + 1)
compute_L <- function(x, s, b){
  one = rep(1,length(b))
  M1 = outer(s^x, b, "*")
  xi_add_one = outer(x,one, "+")
  M2 = outer(s,b, "+")^xi_add_one
  return(M1/M2)
}

## compute ebpm_exponential_mixture problem
ebpm_exponential_mixture <- function(x,s,m = 2, grid = NULL, seed = 123){
  set.seed(seed)
  if(is.null(grid)){grid <- select_grid(x,s,m)}
  b =  grid$b
  L <-  compute_L(x,s,b)
  fit <- mixsqp(L, control = list(verbose = F))
  pi = fit$x
  one = rep(1,length(b))
  cpm = outer(x,one)/outer(s, b)
  lam_pm = cpm %*% pi
  ll = sum(dpois(x, s*lam_pm, log = T))
  return(list(pi = pi, lam_pm = lam_pm, ll = ll,L = L,b = b))
}

## simulate a poisson mean problem
simulate_pm  <-  function(grid = NULL, seed = 123){
  set.seed(seed)
  if(is.null(grid)){lam_true <- c(rep(0,199),seq(0,10,0.05)); pi = NULL}
  else{
    n_component = length(grid$a)
    pi <- rdirichlet(1,rep(1/n_component, n_component))
    lam_true <- rep(sim_mgamma(grid, pi), 400)
  }
  s = rep(1, length(lam_true))
  x  = rpois(length(lam_true),s*lam_true)
  ll = sum(dpois(x, s*lam_true, log = T))
  return(list(x =  x, s = s, lam_true = lam_true, pi = pi, ll = ll))
}

## test functions  above
main <- function(){
  m = 1.1
  sim = simulate_pm()
  x = sim$x 
  s = sim$s
  lam_true = sim$lam_true
  print(sprintf("mle    ll: %f", sum(dpois(x, x, log = T))))
  print(sprintf("oracle ll: %f", sim$ll))
  
  start = proc.time()
  fit = ebpm_exponential_mixture(x, s, m)
  runtime = proc.time() - start
  print(sprintf("fitted ll: %f", fit$ll))
  
  print(sprintf("fit with %d data points and %d grid points", length(sim$x),length(fit$b)))
  print(sprintf("runtime: %f", runtime[[3]]))
  
  df = data.frame(n = 1:length(x), x = x, lam_true = lam_true, lam_hat = fit$lam_pm)
  # ggplot(df)  + geom_point(aes(x = n, y = log10(lam_true+1), color = "true"), cex = 0.5) +
  #   geom_point(aes(x = n, y = log10(x+1), color = "x"), cex = 0.5) +
  #   geom_point(aes(x = n, y = log10(lam_hat+1), color = "lam_hat"), cex = 0.5) +
  #   labs(x = "index", y = "log10(lam + 1)", title = "EBPM") +
  #   guides(fill = "color")
  ggplot(df)  + geom_point(aes(x = x, y = lam_hat, color = "blue"), cex = 0.5) +
    labs(x = "x", y = "lam_hat", title = "EBPM") +
    guides(fill = "color")+
    geom_abline(slope = 1)
}

main()



