#' @title Empirical Bayes Poisson Mean with Point Gamma  as Prior
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Point Gamma: g()  = pi_0 delta() + (1-pi_0) gamma(a, b)
#' 
#' @import stats

#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, a, b)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param g_init vector of initialization for c(pi, a, b); default is c(0.5,1,1)
#' @param control A list of control parameters  to be passed to the optimization function. `nlm` is  used. 
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



# ebpm_point_gamma <- function(x, s, g_init = NULL, seed = 123){
#   if(is.null(g_init)){
#     mus = 10^seq(0,5,1)
#     for(mu in mus){
#       g_init = c(0.5, mu*mean(x[x!=0]/s[x != 0]),mu);
#       try(
#         return(ebpm_point_gamma_helper(x, s, g_init, seed = seed)), silent = T
#       )
#     }
#     stop('ebpm_point_gamma fails  for multiple initialization')
#   }
#   else{
#     try(return(ebpm_point_gamma_helper(x, s, g_init, seed = seed)))
#     stop('ebpm_point_gamma fails  for multiple initialization')
#   }
# }
  



## TODO:
## consider the case where X are all  0s.   If  s are  not all 0, then pi* = 1, which is not reachable after transformation
ebpm_point_gamma <- function(x, s, g_init = NULL, fix_g = F, control = NULL, seed = 123){
  set.seed(seed) ## though seems determined
  #if(!is.null(g_init)){browser()}
  if(is.null(g_init)){g_init = c(0.5,1,1); fix_g =  F}
  g_init = as.numeric(g_init)
  
  if(!fix_g){
    if(is.null(control)){control = nlm_control_defaults()}
    ## MLE
    fn_params = list(x = x, s = s)
    if(all(x  > 0)){ ## in  this case, optimal pi0 is 0, but is not reachable after a transformation in nlm
      g_init = c(g_init[2],  g_init[3])
      opt = do.call(nlm, c(list(pg_nlm_fn_pi0, transform_param_pi0(g_init)), fn_params, control))
      log_likelihood =  -pg_nlm_fn_pi0(opt$estimate, x, s)
      opt_g = c(0, transform_param_back_pi0(opt$estimate))
    }else{
      opt = do.call(nlm, c(list(pg_nlm_fn, transform_param(g_init)), fn_params, control))
      log_likelihood =  -pg_nlm_fn(opt$estimate, x, s)
      opt_g = transform_param_back(opt$estimate)
    }
  }else{
    opt_g = g_init
    log_likelihood =  ifelse(g_init[1] == 0, -pg_nlm_fn_pi0(transform_param_pi0(opt_g), x, s), -pg_nlm_fn(transform_param(opt_g), x, s))
  }
  
  fitted_g = list(pi = opt_g[1], a = opt_g[2], b  = opt_g[3])
  ## posterior mean
  pi = fitted_g$pi
  a =  fitted_g$a
  b =  fitted_g$b
  nb = exp(dnbinom_cts_log_vec(x, a, prob = b/(b+s)))
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
  par0[1] = exp(par[1])
  par0[2] = exp(par[2])
  return(par0)
}

