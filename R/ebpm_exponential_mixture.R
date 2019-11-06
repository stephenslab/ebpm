#' @title Empirical Bayes Poisson Mean with Point Exponential  as Prior
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Mixture of Exponential: \deqn{g()  = sum_k pi_k gamma(shape = 1, rate = b_k)} 
#' b_k is selected to cover the lambda_i  of interest for all data  points x_i
#' @import mixsqp

#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x A vector of Poisson observations.
#' @param s A vector of scaling factors for Poisson observations: the model is \eqn{y[j]~Pois(s[j]*lambda[j])}.
#' @param scale Either \code{"estimate"} if the scale parameters are to be estimated
#'   from the data, or A list of  \code{shape} (set to all 1s for exponential mixture) and  \code{scale} that specifies the  components of the mixture, where each pair of \code{shape}, \code{scale} described a \code{gamma(shape, scale)}. 
#' @param g_init The prior distribution \eqn{g}, of the class \code{gammamix}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. 
#'
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
#'   of estimating it.
#'   
#' @param m multiple coefficient when selectig grid, so the  b_k is of the form {low*m^{k-1}}; must be greater than 1; default is 2
#' @param control A list of control parameters  to be passed to the optimization function. `mixsqp` is  used.
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, and to add posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}}, of class \code{gammamix}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'      }
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    m = 2
#'    out = ebpm::ebpm_exponential_mixture(x,s)
#'    
#' @export

## compute ebpm_exponential_mixture problem
ebpm_exponential_mixture <- function(x,s = 1,  scale = "estimate", g_init = NULL, fix_g = FALSE,m = 2, control =  NULL, low = NULL){
  if(length(s) == 1){s = replicate(length(x),s)}
  if(is.null(control)){control = mixsqp_control_defaults()}
  if(is.null(g_init)){
    fix_g = FALSE ## then automatically unfix g if specified so
    if(identical(scale, "estimate")){scale <- select_grid_exponential(x,s,m, low)}
    g_init = scale2gammamix_init(scale)
  }
  
  if(!fix_g){ ## need to estimate g_hat
    b = 1/g_init$scale ##  from here use gamma(shape = a, rate = b)  where E = a/b
    a = g_init$shape
    tmp <-  compute_L(x,s,a, b)
    L =  tmp$L
    l_rowmax = tmp$l_rowmax
    fit <- mixsqp(L, x0 = g_init$pi, control = control)
    pi = fit$x
    pi = pi/sum(pi) ## seems that some times pi does not sum to one
  }
  else{
    pi = g_init$pi
    a = g_init$shape
    b = 1/g_init$scale
    ## compute loglikelihood
    tmp <-  compute_L(x,s,a, b)
    L =  tmp$L
    l_rowmax = tmp$l_rowmax
  }
  fitted_g = gammamix(pi = pi, shape = a,  scale  = 1/b)
  log_likelihood = sum(log(exp(l_rowmax) * L %*%  pi))
  
  cpm = outer(x,a,  "+")/outer(s, b, "+")
  Pi_tilde = t(t(L) * pi)
  Pi_tilde = Pi_tilde/rowSums(Pi_tilde)
  lam_pm = rowSums(Pi_tilde * cpm)
  c_log_pm = digamma(outer(x,a,  "+")) - log(outer(s, b, "+"))
  lam_log_pm = rowSums(Pi_tilde * c_log_pm)
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(list(fitted_g = fitted_g, posterior = posterior,log_likelihood = log_likelihood))
}

geom_seq <- function(low, up, m){
  N =  ceiling(log(up/low)/log(m)) + 1
  out  = low*m^(seq(1,N, by = 1)-1)
  return(out)
}

# lin_seq <- function(low, up, m){
#   out = seq(low, up, length.out = m)
#   return(out)
# }

## select grid for b_k
select_grid_exponential <- function(x, s, m = 2, d = NULL, low = NULL){
  ## mu_grid: mu =  1/b is the exponential mean
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min(xprime/s)
  mu_grid_max = 2*max(x/s)
  if(is.null(m)){
    if(is.null(d)){m = 2}
    else{m = ceiling((mu_grid_max/mu_grid_min)^(1/(d-1)))}
  }
  if(!is.null(low)){mu_grid_min = min(low, mu_grid_min)}
  mu_grid = geom_seq(mu_grid_min, mu_grid_max, m)
  a = rep(1, length(mu_grid))
  return(list(shape = a, scale = mu_grid))
}



#' @export get_uniform_mixture
get_uniform_mixture <- function(x, s, grid_res = NULL, m = 2, low = NULL){
  if(is.null(grid_res)){
    grid_res = select_grid_exponential(x = x, s = s, m = m, low = low)
  }
  shape = grid_res$shape
  scale = grid_res$scale
  n = length(shape)
  pi = replicate(n, 1/n)
  g = gammamix(pi = pi, shape = shape, scale = scale)
  return(g)
}

## compute L matrix from data and selected grid
## L_ik = NB(x_i; a_k, b_k/b_k + s_i)
## but for computation in mixsqr, we can simplyfy it for numerical stability
#' @export compute_L
compute_L <- function(x, s, a, b){
  prob = 1 - s/outer(s,b, "+")
  l = dnbinom_cts_log(x,a,prob = prob) ## 
  l_rowmax  = apply(l,1,max)
  L = exp(l -  l_rowmax)
  return(list(L = L, l_rowmax = l_rowmax))
}

# it is equivalent to dnbinom in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dnbinom_cts_log <- function(x, a, prob){
  tmp = x*log(1-prob) 
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  out = t(t(log(prob)) * a) + tmp + lgamma(outer(x, a, "+")) - lgamma(x+1)
  out = t(t(out) - lgamma(a))
  return(out)
}

dnbinom_cts_log_1d <- function(x, a, prob){
  tmp = x*log(1-prob)
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  return(a*log(prob) + tmp + lgamma(x+a) - lgamma(x+1) - lgamma(a))
}


