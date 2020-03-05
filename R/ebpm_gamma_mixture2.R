#' @title Empirical Bayes Poisson Mean with Mixture of Gamma as Prior (still in development)
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Mixture of Gamma: \deqn{g()  = sum_k pi_k gamma(shape = a_k, rate = b_k)} 
#' @import mixsqp

#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x A vector of Poisson observations.
#' @param s A vector of scaling factors for Poisson observations: the model is \eqn{y[j]~Pois(s[j]*lambda[j])}.
#' @param shape A vector specifying the shapes used in gamma mixtures
#' @param scale A vector specifying the scales used in gamma mixtures
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
#'         means, and posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}}, of class \code{gammamix}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'      }
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    m = 2
#'    out = ebpm::ebpm_gamma_mixture(x,s)
#'    
#' @export

## compute ebpm_gamma_mixture problem
ebpm_gamma_mixture2 <- function(x,s, grid = NULL,  g_init = NULL, fix_g = FALSE, control_select_grid = NULL, control =  NULL){
  ## a quick  fix when all `x` are 0
  if(max(x) == 0){
    stop("all x are 0") ## TODO replace with sth like gamma(0, ..)
  }
  if(length(s) == 1){s = replicate(length(x),s)}
  if(is.null(control)){control = mixsqp_control_defaults()}
  if(is.null(control_select_grid)){
    ## TODO
  }
  if(is.null(g_init)){
    fix_g = FALSE ## then automatically unfix g if specified so
    if(is.null(grid)){
      params_ = c(list(x = x, s = s), control_select_grid)
      grid = do.call(select_grid2, params_)
    }
    g_init = grid2gammamix(grid, pi = NULL)
  }
  grid = gammamix2grid(g_init) ## make sure g_init and grid are consistent
  tmp <-  compute_L_from_grid(x,s,grid)
  L =  tmp$L
  l_rowmax = tmp$l_rowmax

  ## compute weight pi
  if(!fix_g){ ## need to estimate g_hat
    fit <- mixsqp(L, x0 = g_init$pi, control = control)
    w = fit$x
    w = w/sum(w) ## seems that some times pi does not sum to one
  }
  else{w = g_init$pi}

  fitted_g = grid2gammamix(grid, w)
  log_likelihood = sum(log(exp(l_rowmax) * L %*%  w))
  posterior = compute.posterior.gammamix(x,s,fitted_g, L)
  return(list(fitted_g = fitted_g, posterior = posterior,log_likelihood = log_likelihood))
}

select_grid2 <- function(x, s, mus = NULL , vars = NULL, k = 2){
  if(is.null(mus)){
    if(is.null(k)){ stop("need to provide k in select_grid2") }
    mus = as.vector(kmeans(x, centers = k, nstart = 100, iter.max = 100)$centers)
  }
  if(is.null(vars)){vars = 10^seq(-5,5,1)}
  grid = construct_grid(mus, vars)
  return(grid)
}

construct_grid <- function(mus, vars){
  M = length(mus)
  D = length(vars)
  a = c()
  b = c()
  for(m in 1:M){
    for(d in 1:D){
      b_ = mus[m]/vars[d]
      a = c(a, b_ * mus[m])
      b = c(b, b_)
    }
  }
  return(list(a = a, b = b))
}

compute_L_from_grid <- function(x,s,grid){
  ## TODO: need to consider numerical issue later
  return( compute_L(x = x,s = s,a = grid$a, b = grid$b) )
}

grid2gammamix <- function(grid, pi = NULL){
  n = length(grid$a)
  if(is.null(pi)){pi = replicate(n, 1/n)}
  return( gammamix(pi = pi, shape = grid$a, scale = 1/grid$b) ) 
}

gammamix2grid <- function(g){
  return(list(a = g$shape, b = 1/g$scale))
}

#' @export compute.posterior.gammamix
compute.posterior.gammamix <- function(x,s,g, L){
  a = g$shape
  b = 1/g$scale
  if(is.null(L)){L = compute_L(x = x,s = s,a = g$shape, b = 1/g$scale)$L}
  cpm = outer(x,a,  "+")/outer(s, b, "+")
  Pi_tilde = t(t(L) * g$pi)
  Pi_tilde = Pi_tilde/rowSums(Pi_tilde)
  lam_pm = rowSums(Pi_tilde * cpm)
  c_log_pm = digamma(outer(x,a,  "+")) - log(outer(s, b, "+"))
  lam_log_pm = rowSums(Pi_tilde * c_log_pm)
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(posterior)
}