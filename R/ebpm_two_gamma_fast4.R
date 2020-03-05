# #' @title Empirical Bayes Poisson Mean (mixture of two gammas as prior, faster version 4)
# #' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
# #' with Point Gamma: g()  = pi_0 gamma(shape1, scale1) + (1-pi_0) gamma(shape2, scale2)
# #' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, shape, scale)
# #' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
# #' @param x vector of Poisson observations.
# #' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
# #' @param g_init The prior distribution \eqn{g}, of the class \code{two_gamma}. Usually this is left
# #'   unspecified (\code{NULL}) and estimated from the data. However, it can be
# #'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
# #'   example, to do computations with the "true" \eqn{g} in simulations). If
# #'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
# #'   specifies the initial value of \eqn{g} used during optimization. 
# #' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
# #'   of estimating it.
# #' @param n_iter: number of maximum EM steps
# #' @param rel_tol: tolerance for (maximum) relative change in parameters in \eqn{g}
# #' @param control A list of control parameters to be passed to `mle_two_gamma` function. 
# #' @param get_progress record log_likelihood if set to TRUE
# #' @param seed seed set to initialization if \code{g_init = NULL}

# #' @return A list containing elements:
# #'     \describe{
# #'       \item{\code{posterior}}{A data frame of summary results (posterior
# #'         means, posterior log mean).}
# #'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}} of class \code{point_gamma}} 
# #'       \item{\code{log_likelihood}}{The optimal log likelihood attained
# #'         \eqn{L(\hat{g})}.}
# #'      }

# #' @export
# ebpm_two_gamma_fast4 <- function(x, s = 1, g_init = NULL, 
# 																fix_g = FALSE, n_iter = 100, 
# 																verbose = FALSE, control = NULL, 
# 																get_progress = TRUE, seed = 123){
# 	n = length(x)
#   #if(length(s) == 1){s = replicate(n,s)}
#   # if(is.null(control)){control = tg_control_defaults()}
#   if(is.null(control)){
# 		control = list(nlm_setting = nlm_control_defaults(),
# 								gradient = TRUE, hessian = FALSE)
# 	}
#   if(is.null(g_init)){g_init <- init_two_gamma(x, s, seed = seed);fix_g =  FALSE}
  
#   if(fix_g){
#   	fitted_g = g_init
#   	progress = NULL
#   	log_likelihood = loglik.tg(x, s, fitted_g)
#   	}else{
#   		fit <- mle_two_gamma(x = x, s = s, g = g_init, 
#   	maxiter = n_iter, control = control, verbose = verbose, get_progress = get_progress)
#   		fitted_g = fit$fitted_g
#   		progress = fit$progress
#   		log_likelihood = progress[length(progress)]
#   }
#   ## compute posterior for lambda
#   posterior <- lam.posterior(x, s, fitted_g)
#   return(list(fitted_g = fitted_g, posterior = posterior, log_likelihood = log_likelihood, progress = progress))

# }

# tg_control_defaults <- function(){
# 	list(tol_in = 1e-4)
# }

# init_two_gamma <- function(x, s, seed = 123){
# 	set.seed(seed)
#   #browser()
#   ## use k-means to find 2 clusters
#   ## try up to 5 times 
#   max_try = 5
#   i = 0
#   clst_class = "try-error"
#   while(clst_class == "try-error" & i < max_try){
#   	clst = try(kmeans(x = x/s, centers = 2))
#   	clst_class = class(clst)
#   	i <- i + 1
#   }
  
#   if(class(clst) == "try-error"){ ## then probably there should be only 1 cluster
#   	warning(sprintf("cannot find cluster using k-means after trying %d times", max_try))
#     pi0 = 0
#     shape1 = 1; scale1 = 1;
#     shape2 = 1; scale2 = 1;
#   }else{
#     ## initialzie pi0
#     pi0 = sum(clst$cluster == 1)/length(x)
#     ## estimate shape1, scale1
#     idx = which(clst$cluster == 1)
#     s_ = ifelse(length(s) == 1, s, s[idx])
#     fit_ = ebpm_point_gamma(x = x[idx], s = s_, pi0 = 0)
#     shape1 = fit_$fitted_g$shape
#     scale1 = fit_$fitted_g$scale
#     ## estimate shape2, scale2
#     idx = which(clst$cluster == 2)
#     s_ = ifelse(length(s) == 1, s, s[idx])
#     fit_ = ebpm_point_gamma(x = x[idx], s = s_, pi0 = 0)
#     shape2 = fit_$fitted_g$shape
#     scale2 = fit_$fitted_g$scale
#   }
#   return(list(pi0 = pi0, shape1 = shape1, scale1 = scale1, shape2 = shape2, scale2 = scale2))
# }


# loglik.tg <- function(x, s, g){
# 	pi1 = g$pi0
# 	pi2 = 1 - pi1
#   a1 = g$shape1
#   b1 = 1/g$scale1
#   a2 = g$shape2
#   b2 = 1/g$scale2
#   p1 = b1/(s + b1)
#   p2 = b2/(s + b2)
#   loglik = sum(log(pi1 * exp(loglik.nb(x, a1, p1)) + pi2 * exp(loglik.nb(x, a2, p2))))
#   return(loglik)
# }

# lam.posterior <- function(x, s, g){
#   n = length(x)
#   pi0 = g$pi0 
#   a1 = g$shape1
#   a2 = g$shape2
#   b1 = 1/g$scale1
#   b2 = 1/g$scale2

#   p1 = 1/(1 + g$scale1 * s)
#   p2 = 1/(1 + g$scale2 * s)
  
#   pi1_ = pi0 * exp(loglik.nb(x, a1, p1))
#   pi1 = pi1_/(pi1_ + (1-pi0) * exp(loglik.nb(x, a2, p2)))
#   pi2 = 1 - pi1
#   lam_pm = pi1 * (a1 + x)/(b1 + s) +  pi2 * (a2 + x)/(b2 + s)
#   ## special case: pi2 = 0, the next part is -Inf ...
#   # lam_log_pm = pi1 *(Rfast::Digamma(a1 + x) - log(b1 + s) ) + pi2 *(Rfast::Digamma(a2 + x) - log(b2 + s) )
#   lam_log_pm1 = replicate(n, 0)
#   lam_log_pm2 = replicate(n, 0)
#   if(length(s) == 1){
#   	nz_mask <- (pi1 != 0)
#   	lam_log_pm1[nz_mask] = pi1[nz_mask] * (Rfast::Digamma(a1 + x[nz_mask]) - log(b1 + s)) 
#   	nz_mask <- (pi2 != 0)
#   	lam_log_pm2[nz_mask] = pi2[nz_mask] * (Rfast::Digamma(a2 + x[nz_mask]) - log(b2 + s)) 
# 	  lam_log_pm = lam_log_pm1 + lam_log_pm2
#   }else{
#   	nz_mask <- (pi1 != 0)
#   	lam_log_pm1[nz_mask] = pi1[nz_mask] * (Rfast::Digamma(a1 + x[nz_mask]) - log(b1 + s[nz_mask])) 
#   	nz_mask <- (pi2 != 0)
#   	lam_log_pm2[nz_mask] = pi2[nz_mask] * (Rfast::Digamma(a2 + x[nz_mask]) - log(b2 + s[nz_mask])) 
# 	  lam_log_pm = lam_log_pm1 + lam_log_pm2
#   }
#   posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
#   return(posterior)
# }