#' @title Empirical Bayes Poisson Mean (mixture of two gammas as prior)
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Point Gamma: g()  = pi_0 gamma(shape1, scale1) + (1-pi_0) gamma(shape2, scale2)
#' 
#' @import stats

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
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    out = ebpm_two_gamma(x,s)
#' @export

ebpm_two_gamma <- function(x, s = 1, g_init = NULL, fix_g = F, n_iter = 100, rel_tol = 1e-10,control = NULL){
  n = length(x)
  if(length(s) == 1){s = replicate(n,s)}
  if(is.null(control)){control = nlm_control_defaults()}
  if(is.null(g_init)){g_init = init_two_gamma(x, s); fix_g =  F}
  
  # fitted = ebpm_two_gamma_util(x, s, n_iter, g_init, fix_g)
  progress = replicate(n_iter, -1e+20)
  
  pi0 = g_init$pi0
  a1 = g_init$shape1
  b1 = 1/g_init$scale1
  a2 = g_init$shape2
  b2 = 1/g_init$scale2
  
  if(fix_g){
    log_likelihood = compute_ll_two_gamma(x, s, pi0, a1, b1, a2, b2)
    param_curr = c(pi0, a1, b1, a2, b2)
    fitted_g = two_gamma(pi0 = pi0, shape1 = a1, scale1 = 1/b1, 
                         shape2 = a2, scale2 = 1/b2)
  }else{
    ## EM updates
    for(i in 1:n_iter){
      ## record previous parameters
      param_prev = c(pi0, a1, b1, a2, b2)
      ### E-step: compute Z | X, pi^0
      w1 = compute_posterior_w(x, s, pi0, a1,b1, a2, b2)
      w2 = 1 - w1
      ### M-step:
      ## update pi0
      pi0 = sum(w1)/n
      ## update a1, b1
      tmp_ab = update_ab(w1, x, s, a1, b1, control)
      a1 = tmp_ab$a
      b1 = tmp_ab$b
      ## update a2, b2
      tmp_ab = update_ab(w2, x, s, a2, b2, control)
      a2 = tmp_ab$a
      b2 = tmp_ab$b
      
      ## record progress
      log_likelihood = compute_ll_two_gamma(x, s, pi0, a1, b1, a2, b2)
      progress[i] = log_likelihood
      param_curr = c(pi0, a1, b1, a2, b2)
      rel_diff = compute_rel_diff(param_prev, param_curr)
      
      if(rel_diff < rel_tol){
        #print(sprintf("converges after iter %d", i))
        progress = progress[1:i]
        break
      }
    }
  }
  
  
  fitted_g = two_gamma(pi0 = pi0, shape1 = a1, scale1 = 1/b1, 
                       shape2 = a2, scale2 = 1/b2)
  ## compute posterior
  posterior = compute_posterior_lam(x, s, fitted_g)
  
  return(list(fitted_g = fitted_g, posterior = posterior, log_likelihood = log_likelihood, progress = progress))
}


#######################################################################
######################### Helper Functions ############################
#######################################################################

## compute NB(x, size = a, prob = p)  
compute_nb <- function(x, a, p)
  exp(dnbinom_cts_log_vec(x, a, p))

compute_weighted_nb_log <- function(w, x, a, p)
  sum( w * dnbinom_cts_log_vec(x,a,p))

compute_ll_two_gamma <- function(x, s, pi1, a1, b1, a2, b2){
  nb1 = compute_nb(x, a = a1, p = b1/(b1 + s))
  nb2 = compute_nb(x, a = a2, p = b2/(b2 + s))
  return(sum(log(pi1*nb1 + (1 - pi1)*nb2)))
}


## compute posterior for w: P(Z | X, pi^0)
compute_posterior_w <- function(x, s, pi1, a1,b1, a2, b2){
  n = length(x)
  ## compute posterior Z | X, pi1^0
  w1 =  pi1 * compute_nb(x, a1, b1/(b1+s))   ## P(Z = 1 | X, pi1^0), not scaled yet
  w2 =  (1 - pi1) * compute_nb(x, a2, b2/(b2+s))   ## P(Z = 1 | X, pi1^0), not scaled yet
  w1 = w1/(w1 + w2)
  return(w1)
}


## update a, b in weighted NB
## max_{a,b} sum_i w_i log NB(x_i, a, b/b+s_i)
update_ab <- function(w, x, s, a, b, control){
  fn_params = list(x = x, s = s,  w = w)
  init_t = c(log(a), log(b))
  opt = do.call(nlm, c(list(obj_w_nb, init_t), fn_params, control))
  log_likelihood =  -obj_w_nb(opt$estimate, x, s, w)
  a = exp(opt$estimate[1])
  b = exp(opt$estimate[2])
  return(list(a = a, b = b))
}

## obj for nlm
## par = c(log(a), log(b))
obj_w_nb <- function(par, x, s, w){
  n = length(x)
  a = exp(par[1])
  b = exp(par[2])
  return(- compute_weighted_nb_log(w, x,a, b/(b + s)))
}

## compute relative differences between v_curr and v_prev
compute_rel_diff <- function(v_prev, v_curr){
  return(sum(abs(v_curr - v_prev))/sum(abs(v_curr + v_prev)))
}


compute_posterior_lam <- function(x,s,fitted_g){
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
}

init_two_gamma <- function(x, s){
  #browser()
  ## use k-means to find 2 clusters
  clst = try(kmeans(x = x/s, centers = 2))
  
  if(class(clst) == "try-error"){ ## then probably there should be only 1 cluster
    pi0 = 0
    shape1 = 1; scale1 = 1;
    shape2 = 1; scale2 = 1;
  }else{
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
  }
  return(list(pi0 = pi0, shape1 = shape1, scale1 = scale1, shape2 = shape2, scale2 = scale2))
}












