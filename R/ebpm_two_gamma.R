#' @title Empirical Bayes Poisson Mean (mixture of two gammas as prior)
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Point Gamma: g()  = pi_0 gamma(shape1, scale1) + (1-pi_0) gamma(shape2, scale2)
#' 
#' @import stats

#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, shape, scale)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param g_init The prior distribution \eqn{g}, of the class \code{point_gamma}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. 
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
#'   of estimating it.
#' @param pi0 Either  \code{"estimate"} which optimizes over pi0 along with  shape and scale, or a number in \code{[0,1]} that fixes pi0
#' @param control A list of control parameters  to be passed to the optimization function. `nlm` is  used. 

#' 
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}} of class \code{point_gamma}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'      }
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    out = ebpm_point_gamma(x,s)
#' @export

ebpm_two_gamma <- function(x, s = 1, g_init = NULL, fix_g = F, pi0 = "estimate",control = NULL){
  if(length(s) == 1){s = replicate(length(x),s)}
  if(is.null(control)){control = nlm_control_defaults()}
  if(is.null(g_init)){g_init = init_two_gamma(x, s); fix_g =  F}
  #browser()
  if(!fix_g){
    ## MLE
    if(identical(pi0, "estimate")){
      if(!all(x  > 0)){
        fn_params = list(x = x, s = s)
        opt = do.call(nlm, c(list(tg_nlm_fn, tg_transform_param(g_init)), fn_params, control))
        log_likelihood =  -tg_nlm_fn(opt$estimate, x, s)
        opt_g = tg_transform_param_back(opt$estimate)
      }else{  ## in  this case, optimal pi0 is 0, but is not reachable after a transformation in nlm; so fix pi0 at 0
        pi0 = 0
        fn_params = list(x = x, s = s, pi0 = pi0)
        opt = do.call(nlm, c(list(tg_nlm_fn_fix_pi0, tg_transform_param_fix_pi0(g_init)), fn_params, control))
        log_likelihood =  -tg_nlm_fn_fix_pi0(opt$estimate, x, s, pi0)
        opt_g = c(pi0, tg_transform_param_back_fix_pi0(opt$estimate))
      }
    }else{ ## fix pi0
      pi0 = as.numeric(pi0)
      fn_params = list(x = x, s = s, pi0 = pi0)
      opt = do.call(nlm, c(list(tg_nlm_fn_fix_pi0, tg_transform_param_fix_pi0(g_init)), fn_params, control))
      log_likelihood =  -tg_nlm_fn_fix_pi0(opt$estimate, x, s, pi0)
      opt_g = c(pi0, tg_transform_param_back_fix_pi0(opt$estimate))
    }
  }else{ ## fix_g = T
    opt_g = as.numeric(g_init)
    # log_likelihood =  ifelse( g_init$pi0 == 0, -tg_nlm_fn_pi0(tg_transform_param_pi0(opt_g), x, s), -tg_nlm_fn(tg_transform_param(opt_g), x, s))
    log_likelihood =  -tg_nlm_fn(tg_transform_param(opt_g), x, s)
  }
  fitted_g = list(pi0 = opt_g[1], shape1 = opt_g[2], scale1  = opt_g[3], shape2 = opt_g[4], scale2  = opt_g[5])
  
  ## compute posterior
  pi0 = fitted_g$pi0 
  a1 = fitted_g$shape1
  b1 = 1/fitted_g$scale1
  a2 = fitted_g$shape2
  b2 = 1/fitted_g$scale2
  
  pi1_ = pi0 * exp(dnbinom_cts_log_vec(x, a1, prob = b1/(b1+s)))
  pi2_ = (1-pi0) * exp(dnbinom_cts_log_vec(x, a2, prob = b2/(b2+s)))
  pi1 = pi1_/(pi1_ + pi2_)
  pi2 = pi2_/(pi1_ + pi2_)
  
  lam_pm = pi1 * (a1 + x)/(b1 + s) +  pi2 * (a2 + x)/(b2 + s)
  lam_log_pm = pi1 *( digamma(a1 + x) - log(b1 + s) ) + pi2 *( digamma(a2 + x) - log(b2 + s) )
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(list(fitted_g = fitted_g, posterior = posterior, log_likelihood = log_likelihood))
}

tg_nlm_fn_fix_pi0 <- function(par, x, s, pi0){
  ## d = NB(x, a, b/(b+s))
  ## return - log(pi0 + (1-pi0) d) if x = 0; 
  ## else return - (log(1-pi0) + log(d))
  a1 = exp(par[1])
  b1  =  exp(par[2])
  a2 = exp(par[3])
  b2  =  exp(par[4])
  d1_log <- dnbinom_cts_log_vec(x, a1, b1/(b1+s))
  d2_log <- dnbinom_cts_log_vec(x, a2, b2/(b2+s))
  out = sum( log( pi0*exp(d1_log)  + (1-pi0)*exp(d2_log) ) )
  return(-out)
}

tg_nlm_fn <- function(par, x, s){
  ## d1 = NB(x, a1, b1/(b1+s)); d2 = NB(x, a2, b2/(b2+s))
  ## return - log(pi0 * d1 + (1-pi0) * d2) 
  #browser()
  pi0 = 1/(1+ exp(-par[1]))
  a1 = exp(par[2])
  b1  =  exp(par[3])
  a2 = exp(par[4])
  b2  =  exp(par[5])
  d1_log <- dnbinom_cts_log_vec(x, a1, b1/(b1+s))
  d2_log <- dnbinom_cts_log_vec(x, a2, b2/(b2+s))
  out = sum( log( pi0*exp(d1_log)  + (1-pi0)*exp(d2_log) ) )
  return(-out)
}


## par0: pi0,shape, scale
## want: logit(pi0), log(shape1), log(1/scale1), log(shape2), log(1/scale2)
tg_transform_param <- function(par0){
  par0 = as.numeric(par0)
  par = rep(0,length(par0))
  par[1] = log(par0[1]/(1-par0[1]))
  par[2] = log(par0[2])
  par[3] = - log(par0[3])
  par[4] = log(par0[4])
  par[5] = - log(par0[5])
  return(par)
}

tg_transform_param_back <- function(par){
  par0 = rep(0,length(par))
  par0[1] = 1/(1+ exp(-par[1]))
  par0[2] = exp(par[2])
  par0[3] = exp(- par[3])
  par0[4] = exp(par[4])
  par0[5] = exp(- par[5])
  return(par0)
}

tg_transform_param_fix_pi0 <- function(par0){
  ## only need to optimize over shape, scale
  par0 = c(par0$shape1, par0$scale1, par0$shape2, par0$scale2)
  par = rep(0,length(par0))
  par[1] = log(par0[1])
  par[2] = - log(par0[2])
  par[3] = log(par0[3])
  par[4] = - log(par0[4])
  return(par)
}

tg_transform_param_back_fix_pi0 <- function(par){
  par0 = rep(0,length(par))
  par0[1] = exp(par[1])
  par0[2] = exp(- par[2])
  par0[3] = exp(par[3])
  par0[4] = exp(- par[4])
  return(par0)
}


init_two_gamma <- function(x, s){
  #browser()
  ## use k-means to find 2 clusters
  clst = kmeans(x = x/s, centers = 2)
  ## initialzie pi0
  pi0 = sum(clst$cluster == 1)/length(x)
  ## estimate shape1, scale1
  idx = which(clst$cluster == 1)
  fit_ = ebpm_point_gamma(x = x[idx], s = s[idx], pi0 = 0)
  shape1 = fit_$fitted_g$shape
  scale1 = fit_$fitted_g$scale
  ## estimate shape2, scale2
  idx = which(clst$cluster == 2)
  fit_ = ebpm_point_gamma(x = x[idx], s = s[idx], pi0 = 0)
  shape2 = fit_$fitted_g$shape
  scale2 = fit_$fitted_g$scale
  return(list(pi0 = pi0, shape1 = shape1, scale1 = scale1, shape2 = shape2, scale2 = scale2))
}












