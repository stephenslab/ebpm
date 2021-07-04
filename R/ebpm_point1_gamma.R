#' @title Empirical Bayes Poisson Mean with Point (at 1) Gamma as Prior
#' 
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Point Gamma: g()  = pi_0 delta_1() + (1-pi_0) gamma(shape, scale)
#' 
#' @import stats
#' @import matrixStats
#' 
#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, shape, scale)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param g_init The prior distribution \eqn{g}, of the class \code{point1_gamma}. Usually this is left
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
#'         means, and posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}} of class \code{point1_gamma}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'      }

#' @export


ebpm_point1_gamma <- function(x, s = 1, g_init = NULL, fix_g = FALSE,control = NULL){
  if(length(s) == 1){s = replicate(length(x),s)}
  if(is.null(control)){control = nlm_control_defaults()}
  if(is.null(g_init)){g_init = point1_gamma(0.5,1,1); fix_g =  FALSE}
  ## compute MLE for g
  fit_tmp = mle_point1_gamma(x = x, s = s, g = g_init, fix_g = fix_g, control = control)
  fitted_g = fit_tmp$fitted_g
  log_likelihood = fit_tmp$log_likelihood
  ## compute posterior 
  posterior = posterior_point1_gamma(x = x, s = s, g = fitted_g)
  return(list(fitted_g = fitted_g, posterior = posterior, log_likelihood = log_likelihood))
}

## output: g, loglikelihood
mle_point1_gamma <- function(x, s, g, fix_g, control){
  if(!fix_g){
    fn_params = list(x = x, s = s)
    opt = do.call(nlm, c(list(p1g_nlm_fn, transform_param_p1g(g)), fn_params, control))
    log_likelihood =  -p1g_nlm_fn(opt$estimate, x, s)
    g = transform_param_back_p1g(opt$estimate)
  }else{ 
    log_likelihood = -p1g_nlm_fn(transform_param_p1g(g), x, s)
  }
  return(list(fitted_g = g, log_likelihood = log_likelihood))
}

posterior_point1_gamma <- function(x, s, g){
  pi0 = g$pi0
  shape = g$shape
  scale = g$scale
  ## compute pi0_hat
  nb_log <- dnbinom_cts_log_vec(x = x, a = shape, prob = 1/(1 + scale*s))
  pois_log <- dpois_cts_log_vec(x = x, lam = s)
  lls <- cbind(pois_log + log(pi0), nb_log + log(1-pi0))
  ll = apply(lls, 1, logSumExp)
  pi0_hat = exp(pois_log + log(pi0) - ll)
  ## compute E(lam) & E(log(lam))
  lam_pm = pi0_hat + (1 - pi0_hat) * (shape + x)/(1/scale + s)
  lam_log_pm = (1 - pi0_hat) * (digamma(shape + x) - log(1/scale + s))
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm, pi0_hat = pi0_hat)
  return(posterior)
}


p1g_nlm_fn <- function(par, x, s){
  g = transform_param_back_p1g(par)
  return(-loglik_point1_gamma(x = x, s = s, g = g))
}

loglik_point1_gamma <- function(x, s, g){
  pi0 = g$pi0
  shape = g$shape
  scale = g$scale
  nb_log <- dnbinom_cts_log_vec(x = x, a = shape, prob = 1/(1 + scale*s))
  pois_log <- dpois_cts_log_vec(x = x, lam = s)
  # ll = log(pi0 * exp(pois_log) + (1-pi0) * exp(nb_log))
  lls <- cbind(pois_log + log(pi0), nb_log + log(1-pi0))
  ll = apply(lls, 1, logSumExp)
  return(sum(ll))
}


## par0: pi0,shape, scale
## want: logit(pi0), log(shape), log(scale)
transform_param_p1g <- function(g){
  par = rep(0,3)
  par[1] = log(g$pi0/(1-g$pi0))
  par[2] = log(g$shape)
  par[3] = log(g$scale)
  return(par)
}

transform_param_back_p1g <- function(par){
  point1_gamma(pi0 = 1/(1+ exp(-par[1])), shape = exp(par[2]), scale = exp(par[3]))
}



