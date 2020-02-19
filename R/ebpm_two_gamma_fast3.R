##


ebpm_two_gamma_fast3 <- function(x, s = 1, g_init = NULL, 
	fix_g = FALSE, n_iter = 100, verbose = FALSE, control = NULL, get_progress = TRUE, seed = 123){
	n = length(x)
  #if(length(s) == 1){s = replicate(n,s)}
  if(is.null(control)){control = tg_control_defaults()}
  if(is.null(g_init)){
  	t <- system.time(
  		g_init <- init_two_gamma(x, s, seed = seed)
  	)
  	print("initialization time")
  	print(t)
  	fix_g =  FALSE
  }
  
  if(fix_g){
  	fitted_g = g_init
  	progress = NULL
  	log_likelihood = loglik.tg(x, s, fitted_g)
  	}else{
  		t = system.time(
  		fit <- mle_two_gamma(x = x, s = s, g = g_init, 
  	maxiter = n_iter, tol_in = control$tol_in, verbose = verbose, get_progress = get_progress)
  		)
  		print("mle fitting time")
  		print(t)
  		#browser()
  		fitted_g = fit$fitted_g
  		progress = fit$progress
  		log_likelihood = progress[length(progress)]
  }
  ## compute posterior for lambda
  t = system.time(
  posterior <- lam.posterior(x, s, fitted_g)
  )
  print("posterior computing time")
  print(t)
  return(list(fitted_g = fitted_g, posterior = posterior, log_likelihood = log_likelihood, progress = progress))

}

tg_control_defaults <- function(){
	list(tol_in = 1e-40)
}

init_two_gamma <- function(x, s, seed = 123){
	set.seed(seed)
  #browser()
  ## use k-means to find 2 clusters
  ## try up to 5 times 
  clst_class = "try-error"
  max_try = 5
  i = 0
  while(clst_class == "try-error" & i < max_try){
  	clst = try(kmeans(x = x/s, centers = 2))
  	clst_class = class(clst)
  	i <- i + 1
  }
  
  if(class(clst) == "try-error"){ ## then probably there should be only 1 cluster
    pi0 = 0
    shape1 = 1; scale1 = 1;
    shape2 = 1; scale2 = 1;
  }else{

    ## initialzie pi0
    pi0 = sum(clst$cluster == 1)/length(x)
    ## estimate shape1, scale1
    idx = which(clst$cluster == 1)
    s_ = ifelse(length(s) == 1, s, s[idx])
    fit_ = ebpm_point_gamma(x = x[idx], s = s_, pi0 = 0)
    shape1 = fit_$fitted_g$shape
    scale1 = fit_$fitted_g$scale
    ## estimate shape2, scale2
    idx = which(clst$cluster == 2)
    s_ = ifelse(length(s) == 1, s, s[idx])
    fit_ = ebpm_point_gamma(x = x[idx], s = s_, pi0 = 0)
    shape2 = fit_$fitted_g$shape
    scale2 = fit_$fitted_g$scale
  }
  return(list(pi0 = pi0, shape1 = shape1, scale1 = scale1, shape2 = shape2, scale2 = scale2))
}


loglik.tg <- function(x, s, g){
	pi1 = g$pi0
	pi2 = 1 - pi1
  a1 = g$shape1
  b1 = 1/g$scale1
  a2 = g$shape2
  b2 = 1/g$scale2
  p1 = b1/(s + b1)
  p2 = b2/(s + b2)
  loglik = sum(log(pi1 * exp(loglik.nb(x, a1, p1)) + pi2 * exp(loglik.nb(x, a2, p2))))
  return(loglik)
}

lam.posterior <- function(x, s, g){
  n = length(x)
  pi0 = g$pi0 
  a1 = g$shape1
  a2 = g$shape2
  b1 = 1/g$scale1
  b2 = 1/g$scale2

  p1 = 1/(1 + g$scale1 * s)
  p2 = 1/(1 + g$scale2 * s)
  
  pi1_ = pi0 * exp(loglik.nb(x, a1, p1))
  pi1 = pi1_/(pi1_ + (1-pi0) * exp(loglik.nb(x, a2, p2)))
  pi2 = 1 - pi1
  # lam_pm = pi1 * (a1 + x)/(b1 + s) +  pi2 * (a2 + x)/(b2 + s)
  # lam_log_pm = pi1 *(Rfast::Digamma(a1 + x) - log(b1 + s) ) + pi2 *(Rfast::Digamma(a2 + x) - log(b2 + s) )
  lam_pm = pi1 * (a1 + x)/(b1 + s) +  pi2 * (a2 + x)/(b2 + s)
  # lam_log_pm = pi1 *(Rfast::Digamma(a1 + x) - log(b1 + s) ) + pi2 *(Rfast::Digamma(a2 + x) - log(b2 + s) )
  ## special case: pi2 = 0, the next part is -Inf ...
  lam_log_pm1 = replicate(n, 0)
  lam_log_pm2 = replicate(n, 0)
  if(length(s) == 1){
  	nz_mask <- (pi1 != 0)
  	lam_log_pm1[nz_mask] = pi1[nz_mask] * (Rfast::Digamma(a1 + x[nz_mask]) - log(b1 + s)) 
  	nz_mask <- (pi2 != 0)
  	lam_log_pm2[nz_mask] = pi2[nz_mask] * (Rfast::Digamma(a2 + x[nz_mask]) - log(b2 + s)) 
	  lam_log_pm = lam_log_pm1 + lam_log_pm2
  }else{
  	nz_mask <- (pi1 != 0)
  	lam_log_pm1[nz_mask] = pi1[nz_mask] * (Rfast::Digamma(a1 + x[nz_mask]) - log(b1 + s[nz_mask])) 
  	nz_mask <- (pi2 != 0)
  	lam_log_pm2[nz_mask] = pi2[nz_mask] * (Rfast::Digamma(a2 + x[nz_mask]) - log(b2 + s[nz_mask])) 
	  lam_log_pm = lam_log_pm1 + lam_log_pm2
  }
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(posterior)
}