#' @title Empirical Bayes Poisson Mean with Point Gamma  as Prior
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Point Gamma: g()  = pi_0 delta() + (1-pi_0) gamma(a, b)
#' 
#' @import stats

#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, a, b)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param init_par vector of initialization for c(pi, a, b); default is c(0.5,1,1)
#' @param seed set seed (not necessary now)
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
#'    out = ebpm_point_gamma(x,s)
#' @export


#trace("pg_nlm_fn", quote(if(any(is.nan(log(par[1])))) { browser() }), at=4, print=F)



ebpm_point_gamma <- function(x, s, init_par = NULL, seed = 123){
  if(is.null(init_par)){
    mus = 10^seq(0,5,1)
    for(mu in mus){
      init_par = c(0.5, mu*mean(x[x!=0]/s[x != 0]),mu);
      try(
        return(ebpm_point_gamma_helper(x, s, init_par, seed = seed)), silent = T
      )
    }
    stop('ebpm_point_gamma fails  for multiple initialization')
  }
  else{
    try(return(ebpm_point_gamma_helper(x, s, init_par, seed = seed)))
    stop('ebpm_point_gamma fails  for multiple initialization')
  }
}
  

ebpm_point_gamma_helper <- function(x, s, init_par, seed = 123){
  set.seed(seed) ## though seems determined
  # if(is.null(init_par)){init_par = c(0.5,1,1)}
  # init_par[1] = max(init_par[1], 0.999)
  
  ## MLE
  if(all(x  > 0)){ ## in  this case, optimal pi0 is 0, but is not reachable  after a transformation in nlm
    init_par = c(init_par[2],  init_par[3])
    opt = nlm(pg_nlm_fn_pi0, transform_param_pi0(init_par), x, s, print.level = 0, gradtol = 1e-15)
    log_likelihood =  -pg_nlm_fn_pi0(opt$estimate, x, s)
    opt_par = c(0, transform_param_back_pi0(opt$estimate))
    fitted_g = list(pi = opt_par[1], a = opt_par[2], b  = opt_par[3])
  }else{
    opt = nlm(pg_nlm_fn, transform_param(init_par), x, s, print.level = 0, gradtol = 1e-15)
    log_likelihood =  -pg_nlm_fn(opt$estimate, x, s)
    opt_par = transform_param_back(opt$estimate)
    fitted_g = list(pi = opt_par[1], a = opt_par[2], b  = opt_par[3])
  }
  
  ## posterior mean
  pi = opt_par[1]
  a =  opt_par[2]
  b =  opt_par[3]
  nb = exp(dnbinom_cts_log_vec(x, a, prob = b/(b+s)))
  # lam_pm = ((1-pi)*nb*(a+x)/(b+s))/(pi*as.integer(x ==  0) + (1-pi)*nb)
  pi_hat = pi*as.integer(x ==  0)/(pi*as.integer(x ==  0) + (1-pi)*nb)
  lam_pm = (1-pi_hat)*(a+x)/(b+s)
  lam_log_pm =  digamma(a + x) - log(b + s)
  lam_log_pm[x == 0] = -Inf
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(list(fitted_g = fitted_g, posterior = posterior, log_likelihood = log_likelihood))
}

pg_nlm_fn_pi0 <- function(par, x, s){  ## for the case where x  > 0 for all a, and optimal pi0 is 0
  a  =  exp(par[1])
  b  =  exp(par[2])
  d <- exp(dnbinom_cts_log_vec(x, a, b/(b+s)))
  #if(is.nan(sum(log(pi * c + d)))){browser()}
  return(-sum(log(d)))
}


pg_nlm_fn <- function(par, x, s){
  #browser()
  pi = 1/(1+ exp(-par[1]))
  a = exp(par[2])
  b  =  exp(par[3])
  d <- exp(dnbinom_cts_log_vec(x, a, b/(b+s)))
  c = as.integer(x ==  0) - d
  #if(is.nan(sum(log(pi * c + d)))){browser()}
  return(-sum(log(pi*c + d)))
}

# it is equivalent to dnbinom in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dnbinom_cts_log_vec <- function(x, a, prob){
  tmp = x*log(1-prob)
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  return(a*log(prob) + tmp + lgamma(x+a) - lgamma(x+1) - lgamma(a))
}

transform_param <- function(par0){
  par = rep(0,length(par0))
  par[1] = log(par0[1]/(1-par0[1]))
  par[2] = log(par0[2])
  par[3] = log(par0[3])
  return(par)
}

transform_param_back <- function(par){
  par0 = rep(0,length(par))
  #par0[1] = log(par[1]) - log(1-par[1])
  par0[1] = 1/(1+ exp(-par[1]))
  par0[2] = exp(par[2])
  par0[3] = exp(par[3])
  return(par0)
}

transform_param_pi0 <- function(par0){
  par = rep(0,length(par0))
  par[1] = log(par0[1])
  par[2] = log(par0[2])
  return(par)
}

transform_param_back_pi0 <- function(par){
  par0 = rep(0,length(par))
  #par0[1] = log(par[1]) - log(1-par[1])
  par0[1] = exp(par[1])
  par0[2] = exp(par[2])
  return(par0)
}

