#' @title Empirical Bayes Poisson Mean (gh-gamma as prior)
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{s_j lambda_j ~ g()}
#' with g() being gh-gamma
#' @import gsl
#' @import numDeriv
#' @export

ebpm_gh_gamma <- function(x, s = 1, g_init = NULL, fix_g = c(F,F,F,F,F),control = NULL){
  if(is.null(control)){control = nlm_control_defaults()}
  ## initialize g
  tmp_init <- init.gh(x, g_init, fix_g)
  ## get g(theta), theta | x, and log p(x | g)
  fit_theta = ebpm_gh_gamma_util(x = x, var.est = tmp_init$var.est, var.fix = tmp_init$var.fix, control = control)
  ## convert results for theta to that for lam, where theta = s * lam (note there is no clear way for converting g, so I keep them unchanged)
  posterior = list(mean = fit_theta$posterior$mean/s, mean_log = fit_theta$posterior$mean_log - log(s))
  return(list(fitted_g = fit_theta$fitted_g, log_likelihood = fit_theta$log_likelihood, posterior = posterior))
}


## x ~ Po(theta)
## theta ~ Ga(alpha, k/1-k)
## k ~ GH(a,b,phi, gam)
## so g() = g(.;a,b,alpha, phi, gam)
ebpm_gh_gamma_util <- function(x, var.est, var.fix, control){
  var.n = c(names(var.est), names(var.fix))
  opt =  mle.gh.gamma(y = x, var.fix = var.fix, var.est = var.est, control = control)
  
  log_likelihood = - opt$minimum
  
  fitted_g = as.list(c(opt$estimate, var.fix))
  names(fitted_g) = var.n
  fitted_g = list2gh(fitted_g)
  
  posterior = compute_posterior.gh.theta(x, fitted_g)
  
  return(list(fitted_g = fitted_g, log_likelihood = log_likelihood, posterior = posterior))
}

