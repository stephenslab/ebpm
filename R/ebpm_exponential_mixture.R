#' @title Empirical Bayes Poisson Mean with Point Gamma  as Prior
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Mixture of Exponential: g()  = \sum_k \pi_k gamma(1, b_k) 
#' b_k is selected to cover the lambda_i  of interest for all data  points x_i
#' @import mixsqp


#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param m multiple coefficient when selectig grid, so the  b_k is of the form {low*m^{k-1}}; must be greater than 1; default is 2
#' @param grid locations of b_k; if left NULL the algorithm automatically selects them based on data
#' @param seed random seed
#' 
#' 
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, and to add posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'       \item{\code{posterior_sampler}}{(TO ADD!!!) A function that can be used to
#'         produce samples from the posterior. It takes a single parameter
#'         \code{nsamp}, the number of posterior samples to return per
#'         observation.}
#'      }
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    m = 2
#'    out = ebpm_exponential_mixture(x,s, m)
#'    
#' @export

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
  log_likelihood = sum(log(exp(l_rowmax) * L %*%  fit$x))
  pi = fit$x
  fitted_g = list(pi = pi, a = a,  b  = b)
  cpm = outer(x,a,  "+")/outer(s, b, "+")
  Pi_tilde = t(t(L) * pi)
  Pi_tilde = Pi_tilde/rowSums(Pi_tilde)
  lam_pm = rowSums(Pi_tilde * cpm)
  posterior = list(mean = lam_pm)
  return(list(fitted_g = fitted_g, posterior = posterior,log_likelihood = log_likelihood))
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