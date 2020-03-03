## mle for twogamma-poisson model
## x_i \sim Pois(s \lambda_i)
## \lambda_i \sim \pi_0 Ga(shape1, scale1) + (1 - \pi_0) Ga(shape2, scale2) 
#' @import Rfast
#' @import stats
#' @export
mle_two_gamma <- function(x, s, 
													control, 
													g, maxiter = 1, tol_in = 1e-09, 
													verbose = FALSE, get_progress = TRUE){
	pi1 = g$pi0
	pi2 = 1 - pi1
  a1 = g$shape1
  b1 = 1/g$scale1
  a2 = g$shape2
  b2 = 1/g$scale2
  fit = mle_two_gamma_wh(x = x, s = s, pi1 = pi1, pi2 = pi2, 
													a1 = a1, b1 = b1, a2 = a2, b2 = b2,control = control, 
													maxiter = maxiter, tol_in = tol_in, 
													verbose = verbose, get_progress = get_progress)
  g_ = fit$param
  progress = fit$logliks
  fitted_g = two_gamma(pi0 = g_$pi1, 
  	shape1 = g_$a1, scale1 = 1/g_$b1, 
  	shape2 = g_$a2, scale2 = 1/g_$b2)
  return(list(fitted_g = fitted_g, progress = progress))
}

## mle for mixture of two negative binomial
## x ~ pi1 NB(.; a1, p1) + pi2 NB(.; a2, p2)
## where p1 = b1/(b1 + s), p2 = b2/(b2 + s)
mle_two_gamma_wh <- function(x, s, pi1, pi2, a1, b1, a2, b2, 
													control, 
													maxiter = 10, tol_in = 1e-9, 
													verbose = FALSE, get_progress = TRUE){
	n = length(x)
	if(get_progress){logliks <- c()}
	else{logliks <- NULL}
	if(verbose & get_progress){print("iter 			loglik\n")}

	## EM updates
	for(i in 1:maxiter){
		## E-step
		w1 = pi.posterior(x = x, s = s, pi1 = pi1, pi2 = pi2, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
		w2 = 1 - w1
		w1_scaled = w1/sum(w1)
		w2_scaled = w2/sum(w2)
		## M-step
		## update pi0
		pi1 = sum(w1)/n
		pi2 = 1 - pi1
		### update a1, p1
		fit_res = w.negbin.mle(x = x, s = s, w = w1_scaled, expr1 = a1, tol = tol_in, control = control)
		a1 = fit_res$param$a
		b1 = fit_res$param$b


		### update a2, p2
		fit_res = w.negbin.mle(x = x, s = s, w = w2_scaled, expr1 = a2, tol = tol_in, control = control)
		a2 = fit_res$param$a
		b2 = fit_res$param$b
		## check progress
		if(get_progress){logliks <- 
			c(logliks, sum(log(pi1 * exp(loglik.nb(x = x, s = s, a = a1, b = b1)) + 
			pi2 * exp(loglik.nb(x = x, s = s, a = a2, b = b2)))))}
		if(verbose & get_progress){print(sprintf("%d 		%f\n", i, logliks[i]))}
	}

	loglik <- ifelse(get_progress, logliks[i], 
		sum(log(pi1 * exp(loglik.nb(x = x, s = s, a = a1, b = b1)) + 
		pi2 * exp(loglik.nb(x = x, s = s, a = a2, b = b2)))))
	param = list(pi1 = pi1, a1 = a1, b1 = b1, pi2 = pi2, a2 = a2, b2 = b2)
	return(list(param = param, logliks = logliks, loglik = loglik))
}


## mle for weighted negative binomial
## weight `w` sums to 1
## expr1 is the initiaization for size `a`
#' @import Rfast
w.negbin.mle <- function(x, s, w, control, expr1 = NULL, tol = 1e-05) {
	start = proc.time()
	## initialization for a = expr1, and p = a/(a + m)
	n <- length(x)
	m <- sum(w*x)
	if(is.null(expr1)){	
		m2 <- sum(w * x^2) 
		p_ <- 1 - m / (m2 - m^2)
		expr1 <- m / p_ - m ## not sure!!
	}
	r1 <- log(expr1)
	fn_params = list(x = x, s = s, w = w, m = m, 
								gradient = control$gradient, hessian = control$hessian)
	opt = do.call(nlm, c(list(obj.wnb.loga, r1),fn_params, control$nlm_setting))

	## get output
	r2 = opt$estimate[1]
	expr2 <- exp(r2)
  b <- s * expr2 / m
  param <- c( b, expr2, m )
  names(param) <- c("b", "a", "mean")
  param = as.list(param)
  loglik <- sum(loglik.wnb(x,s, w,expr2,b))
  runtime = proc.time() - start
  return(list(iters = opt$iterations, loglik = loglik, param = param, runtime = runtime[[3]]))
}

obj.wnb.loga <- function(r,x,s, w, m = NULL, gradient = TRUE, hessian = FALSE){
	if(is.null(m)){m <- sum(w*x)}
	expr = exp(r)
	# p <- expr / (expr + m) 
	b = s * expr / m
	out <- -sum(loglik.wnb(x,s, w,expr,b))
	if(gradient){
		g <- gradient.wnb.loga(r, x, w, m)
		attr(out, 'gradient') <- g
	}
	if(hessian){
		attr(out, 'hessian') <- hessian.wnb.loga(r, x, w, g, m)
	}
	return(out)
}


## compute gradient of loglik.nb w.r.t r where r = log(a)
gradient.wnb.loga <- function(r, x, w, m = NULL){
	if(is.null(m)){m = sum(w*x)}
	expr = exp(r)
	# g <- sum( w * diff.psigamma(expr, x, 0) ) * expr  + 
	# 	 expr * r - expr * log(expr + m)  
	## for numerical stability when m << expr
	g <- sum( w * diff.psigamma(expr, x, 0) ) * expr  - expr * log1p(m/expr)
	return(- g)
}

## compute hessian of loglik.nb w.r.t r where r = log(a)
hessian.wnb.loga <- function(r, x, w, g = NULL, m = NULL){
	if(is.null(m)){m = sum(w*x)}
	if(is.null(g)) g = gradient.wnb.loga(x,a,w)
	expr = exp(r) 
	h = - g + sum( w * diff.psigamma(expr, x, 1) ) * expr^2  + expr/(1 + expr/m)
	return(- h)
}




## compute E[z_i= 1|x_i,a,p]
pi.posterior <- function(x, s, pi1, pi2, a1, b1, a2, b2){
	n = length(x)
	w1 = 1 / (1 + (pi2/pi1)  * exp( loglik.nb(x, s, a2, b2) - loglik.nb(x, s, a1, b1)) )
	return(w1)
}

## compute `psigamma(x + dx, deriv) - psigamma(x, deriv)`
## x is s scalar, can be large
## dx is a vector, can be very small

## deriv    type
## -1		lgamma
## 0 		digamma
## 1		trigamma
diff.psigamma <- function(x, dx, deriv){
	n = length(dx)
	diff = replicate(n, 0)

	nz_mask <- (dx != 0)
	mask_taylor = (dx < 1e-3) & (x > 1e+4)
	mask_direct = !mask_taylor & nz_mask
	mask_taylor = mask_taylor & nz_mask

	if(any(mask_direct)){
		f0 = psigamma.my(x, deriv = deriv)
		diff[mask_direct] = psigamma.my(x + dx[mask_direct], deriv = deriv) - f0
	}
	if(any(mask_taylor)){ ## large + small results in numerical issue; use taylor expansion 
		f1 = psigamma(x, deriv = deriv + 1)
		f2 = psigamma(x, deriv = deriv + 2)
		diff[mask_taylor] = f1 * dx[mask_taylor]  + 0.5 * f2 * dx[mask_taylor]^2
	}
	return(diff) 
}


psigamma.my <- function(x, deriv){
	if(deriv == -1){return(Rfast::Lgamma(x))}
	if(deriv == 0){return(Rfast::Digamma(x))}
	if(deriv == 1){return(Rfast::Trigamma(x))}
}


## compute loglikelihood for negative binomial (defined using gamma-poisson mixture)
## `a` is size; `b` is defined so that `b = b/(b + s)` where `p` is probability 
## when `x` is non-integer, use the same PDF
## allow `x` and `p` to be vector; `a` must be scalar
# out = diff.psigamma(a, x, -1) - Rfast::Lgamma(x + 1) + a * log(p) +  x * log( 1- p)
## output is vector of the same size as x (and possibly p)
loglik.nb <- function(x, s, a, b){
	logp_ = - log1p(s/b)
	log1_p_ = log(s) - log(b) + logp_
	## special case:
	## x = 0, p = 1, the loglikelihood is finite
	out = diff.psigamma(a, x, -1) - Rfast::Lgamma(x + 1) + a * logp_ 
	nz_mask <- (x!=0)
	out[nz_mask] = out[nz_mask] +  x[nz_mask] * log1_p_
	return(out)
}

## compute weighted loglikelihood for negative binomial
## input / output format is the same as loglik.nb
# out = w * (diff.psigamma(a, x, -1) - Rfast::Lgamma(x + 1) + a * log(p) +  x * log( 1- p))
loglik.wnb <- function(x,s,w,a,b){
	logp_ = - log1p(s/b)
	log1_p_ = log(s) - log(b) + logp_
	## special case
	## consider w_i * x_i * log( 1- p_i)
	## p_i = 0, (x_i = 0 or w_i = 0) would give 0, not NA
	wx = w*x
	nz_mask <- (wx != 0)
	out = w * (diff.psigamma(a, x, -1) - Rfast::Lgamma(x + 1) + a * logp_)
	out[nz_mask] = out[nz_mask] + wx[nz_mask] * log1_p_
	return(out)
}



#### no gradient

# [1] "compare runtime"
# [1] "tg"
#    user  system elapsed 
#  19.131   0.004  19.135 
# [1] "tg_fast"
#    user  system elapsed 
#   2.996   0.000   2.997 
# [1] "compare loglikelihood"
# [1] "tg"
# [1] -1498.943
# [1] "tg_fast"
# [1] -1498.941