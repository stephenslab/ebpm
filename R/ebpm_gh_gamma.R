##

ebpm_gh_gamma <- function(x, s = 1, g_init = NULL, fix_g = c(F,F,F,F,F), pi0 = "estimate",control = NULL){
  ## initialize g
  tmp_init <- init.gh(x, g_init, fix_g)
  ## get g(theta), theta | x, and log p(x | g)
  fit_theta = ebpm_gh_gamma_util(x = x, var.est = tmp_init$var.est, var.fix = tmp_init$var.fix, control = control)
  ## convert results for theta to that for lam, where theta = s * lam (note there is no clear way for converting g, so I keep them unchanged)
  posterior = list(mean = fit_theta$posterior$mean/s. mean_log = fit_theta$posterior$mean_log - log(s))
  return(list(fitted_g = fit_theta$fitted_g, log_likelihood = fit_theta$log_likelihood, posterior = posterior))
}


## x ~ Po(theta)
## theta ~ Ga(alpha, k/1-k)
## k ~ GH(a,b,phi, gam)
## so g() = g(.;a,b,alpha, phi, gam)
ebpm_gh_gamma_util <- function(x, var.est, var.fix, control){
  opt =  mle.gh.gamma(var.fix = var.fix, var.est = var.est, control = control)
  var.gh = update_var(opt = opt, var.est = var.est, var.fix = var.fix)
  log_likelihood = - opt$minimum
  fitted_g = list2gh(var.gh)
  posterior = compute_posterior.gh.theta(x, fitted_g)
  return(list(fitted_g = fitted_g, log_likelihood = log_likelihood, posterior = posterior))
}

update_var <- function(opt, var.est, var.fix){
  var.est.n = names(var.est)
  var.est = as.list(transform.gh.gamma(opt$estimate, back = TRUE))
  names(var.est) = var.est.n
  return(c(var.est, var.fix))
}

list2gh <- function(var){
  ## todo
}

compute_posterior.gh.theta <- function(x, g){
  posterior.theta.mean = (x + g$alpha) * compute_posterior.1minus_psi(x, g, log = FALSE)
  posterior.theta.mean_log = digamma(x + g$alpha) + compute_posterior.1minus_psi(x, g, log = TRUE)
  return(list(mean = posterior.theta.mean, mean_log = posterior.theta.mean_log))
}

compute_posterior.1minus_psi <- function(x, g, log){
  a = g$a + g$alpha
  b = g$b + x
  gam = g$gam
  z = - (1 - g$phi)
  if(!log){ ## compute E[1 - X]
    beta_div = b/(a + b) ## B(a, b+1)/B(a, b) = a/(a + b), from http://mathworld.wolfram.com/BetaFunction.html
    f21_div = f21.my.vector(a = gam, b = a, c = a + b + 1, z = -z)/f21.my.vector(a = gam, b = a, c = a + b, z = -z)
    return(beta_div * f21_div)
  }else{## compute E[log(1 - X)]
    beta_part = digamma(b) - digamma(a + b)
    ## compute f21 part: f21^(0010)/f21^(0000)
    denom = f21.my.vector(a = gam, b = a, c = a + b, z = -z)
    numer = compute.f21.deriv3(a = gam, b = a, c = a + b, z = -z)
    f21_part = numer / denom
    return(beta_part + f21_part)
  }
}

compute.f21.deriv3 <- function(a,b,c,z){
 ## todo
}



mle.gh.gamma <- function(y, var.fix, var.est, control){
  #browser()
  ## need to consider: var.fix is empty (done); var.est is empty (not yet)
  tmp = process.param(var.est, var.fix)
  fn_params = list(y = y, var.fix = tmp$var.fix, var.n = tmp$var.n)
  opt = do.call(nlm, c(list(obj.gh.gamma, tmp$var.est.t), fn_params, control))
  return(opt)
}

transform.gh.gamma <- function(vars, back = FALSE){
  # return(ifelse(back, yes = exp(vars), no = log(vars)))
  if(back){return(exp(vars))}
  else{return(log(vars))}
}

obj.gh.gamma <- function(var.est.t, y, var.fix, var.n){
  #browser()
  if(length(var.est.t) == 6){vars = transform.gh.gamma(var.est.t, back = TRUE)}
  else{vars = c(transform.gh.gamma(var.est.t, back = TRUE), var.fix)}
  #vars = c(transform.gh.gamma(var.est.t, back = TRUE), var.fix)
  names(vars) = var.n
  obj = - loglik.gh.gamma(y = y, gam = vars[["gam"]], phi = vars[["phi"]], 
                          alpha = vars[["alpha"]], a = vars[["a"]], b = vars[["b"]])
  return(obj)
}

## output a list
process.param <- function(var.est, var.fix){
  var.n = c(names(var.est), names(var.fix))
  var.est.t = transform.gh.gamma(as.numeric(var.est))
  var.fix = as.numeric(var.fix)
  return(list(var.est.t = var.est.t, var.fix = var.fix, var.n = var.n))
}


## this speeds up a lot (compared to hypergeo function) !!
f21.my <- function(A, B,C,z){
  if(z < 1 && z > -1){return(f21hyper(A,B,C,z))}
  return(+Inf)
}

f21.my.vector <- Vectorize(f21.my)

loglik.gh.gamma <- function(y, gam, phi, alpha, a, b){
  tmp = lgamma(y + alpha) - lgamma(y + 1) - lgamma(alpha) ## could be a problem when big number adds small number
  tmp = tmp + lbeta(alpha + a, y + b) - lbeta(a, b) ## again, big number add small number ..
  # tmp = tmp + log(hypergeo(A = gam, B = alpha + a, C = y + alpha + a + b, z = 1 - phi)) - log(hypergeo(A = gam, B = a, C = a + b, z = 1 - phi))
  tmp = tmp + log(f21.my.vector(A = gam, B = alpha + a, C = y + alpha + a + b, z = 1 - phi)) - log(f21.my.vector(A = gam, B = a, C = a + b, z = 1 - phi))
  return(sum(Re(tmp)))
}
