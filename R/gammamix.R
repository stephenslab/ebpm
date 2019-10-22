#' @title \code{gammamix} class
#' @description gammamix class that describes \deqn{g()  = sum_k pi_k gamma(shape, scale)} 
#' @param pi A vector of mixture proportions.
#' @param shape A vector of shapes
#' @param scale A vector of scales
#'
#' @export 
gammamix <- function(pi, shape, scale) {
structure(data.frame(pi = pi, shape = shape, scale = scale), class="gammamix")
}

scale2gammamix_init <- function(scale){
  n = length(scale$shape)
  return(gammamix(pi = replicate(n, 1)/n, shape = scale$shape, scale =  scale$scale))
}