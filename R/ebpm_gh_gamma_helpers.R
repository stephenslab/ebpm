#####################################################################
############################INITIALIZATION###########################
#####################################################################
## input example
## fix_g = list(TRUE,TRUE,TRUE, FALSE,FALSE) or simply fix_g = TRUE; note the order corresponds to g_init!!
## g_init = list(a = 0.5, b = 0.5, alpha = 0.5, phi = NULL, gam = NULL), or of the class `gh_gamma`
## output: list(var.est, var.fix)
init.gh <- function(x, g_init, fix_g){
  if(is.null(g_init)){g_init =  list(a = NULL, b = NULL, alpha = NULL, phi = NULL, gam = NULL)}
  var.n = names(g_init)
  ## fill NULL with default settings
  for(n_ in var.n){
    if(is.null(g_init[[n_]])){
      if(n_ == "a"){g_init[[n_]] = 0.5}
      if(n_ == "b"){g_init[[n_]] = 0.5}
      if(n_ == "alpha"){g_init[[n_]] = 0.5}
      if(n_ == "phi"){g_init[[n_]] = min(0.9, sum(x > 0)/sum(x == 0))} ## need to improve
      if(n_ == "gam"){
        gam = try(mean(kmeans(x = x, centers = 2)$centers))
        if(gam == "try-error"){gam = mean(x)}
        g_init[[n_]] = gam
      }
    }
  }
  var.est = g_init[!fix_g]
  var.fix = g_init[fix_g]
  return(list(var.est = var.est, var.fix = var.fix))
}
#####################################################################
############################ MLE ####################################
#####################################################################
mle.gh.gamma <- function(y, var.fix, var.est, control){
  ## if all vars are fixed 
  if(length(var.est) == 0){
    vars = var.fix
    minimum = - loglik.gh.gamma(y = y, gam = vars[["gam"]], phi = vars[["phi"]],
                                alpha = vars[["alpha"]], a = vars[["a"]], b = vars[["b"]])
    estimate = as.numeric(var.fix)
    return(list(minimum = minimum, estimate = estimate))
  }
  #browser()
  var.n = c(names(var.est), names(var.fix))
  fn_params = list(y = y, var.fix = as.numeric(var.fix), var.n = var.n)
  opt = do.call(nlm, c(list(obj.gh.gamma, as.numeric(var.est)), fn_params, control))
  return(opt)
}

obj.gh.gamma <- function(y, var.est,var.fix, var.n){
  vars = as.list(c(var.est, var.fix))
  names(vars) = var.n
  # inf.my = .Machine$double.xmax
  inf.my = 1e+30
  if(vars$a < 0 || vars$a == 0){return(inf.my)}
  if(vars$b < 0 || vars$b == 0){return(inf.my)}
  if(vars$alpha < 0 || vars$alpha == 0){return(inf.my)}
  if(vars$gam < 0){return(inf.my)} ## TODO: check if gam == 0 will cause problems
  if(vars$phi < 0 || vars$phi ==0 || vars$phi > 2 || vars$phi == 2){return(inf.my)}
  obj = - loglik.gh.gamma(y = y, gam = vars[["gam"]], phi = vars[["phi"]],
                          alpha = vars[["alpha"]], a = vars[["a"]], b = vars[["b"]])
  return(obj)
}
#####################################################################
############################ compute posterior ######################
#####################################################################
compute_posterior.gh.theta <- function(x, g){
  posterior.theta.mean = (x + g$alpha) * compute_posterior.gh.1minus_psi(x, g, log = FALSE)
  posterior.theta.mean_log = digamma(x + g$alpha) + compute_posterior.gh.1minus_psi(x, g, log = TRUE)
  return(list(mean = posterior.theta.mean, mean_log = posterior.theta.mean_log))
}

compute_posterior.gh.1minus_psi <- function(x, g, log){
  a = g$a + g$alpha
  b = g$b + x
  gam = g$gam
  z = - (1 - g$phi)
  if(!log){ ## compute E[1 - X]
    beta_div = b/(a + b) ## B(a, b+1)/B(a, b) = a/(a + b), from http://mathworld.wolfram.com/BetaFunction.html
    f21_div = Gauss2F1(A = gam, B = a, C = a + b + 1, z = -z)/Gauss2F1(A = gam, B = a, C = a + b, z = -z)
    return(beta_div * f21_div)
  }else{## compute E[log(1 - X)]
    beta_div_log = digamma(b) - digamma(a + b)
    ## compute f21_div_log: f21^(0010)/f21^(0000)
    f21_div_log = Gauss2F1.deriv3(A = gam, B = a, C = a + b, z = -z) / Gauss2F1(A = gam, B = a, C = a + b, z = -z)
    return(beta_div_log + f21_div_log)
  }
}
#####################################################################
############################ compute likelihoods ####################
#####################################################################
loglik.gh.gamma <- function(y, gam, phi, alpha, a, b){
  # tmp = lgamma(y + alpha) - lgamma(alpha) - lgamma(y + 1)  ## could be a problem when big number adds small number
  ## compute lgamma(alpha + y) - lgamma(alpha)
  if(alpha > 1e+4 && y < 1e+4 && alpha/y > 1e+4){ ## this condition needs refinement!!
    lgamma_diff = lgamma_diff_taylor(alpha, y)
  }
  else{lgamma_diff = lgamma(y + alpha) - lgamma(alpha) }
  
  tmp = lgamma_diff - lgamma(y + 1)
  tmp = tmp + lbeta(alpha + a, y + b) - lbeta(a, b) ## again, big number add small number ..
  tmp = tmp 
  + log(Gauss2F1(A = gam, B = alpha + a, C = y + alpha + a + b, z = 1 - phi))
  - log(Gauss2F1(A = gam, B = a, C = a + b, z = 1 - phi))
  return(sum(tmp))
}

Gauss2F1 <- function(A,B,C,z){
  if(z>=0 & z<1){
    hyperg_2F1(A,B,C,z)
  }else{
    hyperg_2F1(C-A,B,C,1-1/(1-z))/(1-z)^B
  }
}

Gauss2F1.deriv3 <- function(A,B,C,z){
  # browser()
  return(grad(Gauss2F1, x = C, A = A, B = B, z = z, method = "Richardson"))
}


#####################################################################
############################ others ####################
#####################################################################
list2gh <- function(var){
  gh_gamma(alpha = var$alpha, a = var$a, b = var$b, phi = var$phi, gam = var$gam) 
}
