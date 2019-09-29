#' @title Empirical Bayes Poisson Mean with Point Gamma  as Prior
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Point Gamma: g()  = pi_0 \delta() + (1-\pi_0) gamma(a, b)
#' #' @import stats

#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, a, b)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param init_par vector of initialization for c(pi, a, b); default is c(0.5,1,1)
#' @param seed set seed (not necessary now)
#' 
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    out = ebpm_point_gamma(x,s)
#' @export


ebpm_point_gamma <- function(x, s, init_par = c(0.5,1,1), seed = 123){
  set.seed(seed) ## though seems determined
  ## MLE
  opt = nlm(pg_nlm_fn, transform_param(init_par), x, s)
  opt_par = transform_param_back(opt$estimate)
  ll =  -pg_nlm_fn(transform_param(opt_par), x, s)
  ## posterior mean
  pi = opt_par[1]
  a =  opt_par[2]
  b =  opt_par[3]
  nb = dnbinom(x, size = a, prob = b/(b+s))
  pm = ((1-pi)*nb*(a+x)/(b+s))/(pi*as.integer(x ==  0) + (1-pi)*nb)
  return(list(param = opt_par, lam_pm = pm, ll = ll))
}

pg_nlm_fn <- function(par, x, s){
  pi = log(par[1]) - log(1-par[1])
  a = par[2]
  b  =  exp(par[3])
  d <- dnbinom(x, a, b/(b+s), log = F) 
  c = as.integer(x ==  0) - d
  return(-sum(log(pi*c + d)))
}

transform_param <- function(par0){
  par = rep(0,length(par0))
  par[1] = 1/(1+exp(-par0[1]))
  par[2] = par0[2]
  par[3] = log(par0[3])
  return(par)
}

transform_param_back <- function(par){
  par0 = rep(0,length(par))
  par0[1] = log(par[1]) - log(1-par[1])
  par0[2] = par[2]
  par0[3] = exp(par[3])
  return(par0)
}