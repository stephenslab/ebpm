#' @title Empirical Bayes Poisson Mean (invBeta-gamma as prior)
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{s_j lambda_j ~ g()}
#' with g() being invBeta-gamma
#' @import gsl
#' @import numDeriv
#' @export

#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, shape, scale)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param g_init The prior distribution \eqn{g}, of the class \code{two_gamma}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. 
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
#'   of estimating it.
#' @param n_iter: number of maximum EM steps
#' @param rel_tol: tolerance for (maximum) relative change in parameters in \eqn{g}
#' @param control A list of control parameters  to be passed to the optimization function. `nlm` is  used. 

#' @return A list containing elements:
#'     \describe{
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}} of class \code{point_gamma}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'      }


ebpm_invBeta_gamma <- function(x, s = 1, g_init = NULL, fix_g = c(F,F,F),control = NULL){
  #browser()
	if(is.null(control)){control = nlm_control_defaults()}
  ## initialize g (expand to be like gh-gamma)
  tmp_init <- init.invBeta(x, g_init, fix_g)
  ## get g(theta), theta | x, and log p(x | g)
  fit_theta = ebpm_gh_gamma_util(x = x, var.est = tmp_init$var.est, var.fix = tmp_init$var.fix, control = control)
  ## convert results for theta to that for lam, where theta = s * lam (note there is no clear way for converting g, so I keep them unchanged)
  posterior = list(mean = fit_theta$posterior$mean/s, mean_log = fit_theta$posterior$mean_log - log(s))
  return(list(fitted_g = fit_theta$fitted_g[c("a", "b", "alpha")], log_likelihood = fit_theta$log_likelihood, posterior = posterior))
}

## input example
## fix_g = list(TRUE,TRUE,TRUE, FALSE,FALSE) or simply fix_g = TRUE; note the order corresponds to g_init!!
## g_init = list(a = 0.5, b = 0.5, alpha = 0.5, phi = NULL, gam = NULL), or of the class `gh_gamma`
## output: list(var.est, var.fix)
init.invBeta <- function(x, g_init, fix_g){
  if(is.null(g_init)){g_init =  list(a = NULL, b = NULL, alpha = NULL)}
  var.n = names(g_init)
  ## fill NULL with default settings
  for(n_ in var.n){
    if(is.null(g_init[[n_]])){
      if(n_ == "a"){g_init[[n_]] = 0.5}
      if(n_ == "b"){g_init[[n_]] = 0.5}
      if(n_ == "alpha"){g_init[[n_]] = 0.5}
    }
  }
  var.est = g_init[!fix_g]
  var.fix = g_init[fix_g]
  var.fix = c(var.fix, list(phi = 1, gam = 0))
  return(list(var.est = var.est, var.fix = var.fix))
}