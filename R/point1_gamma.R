#' @title \code{point1_gamma} class
#' @description gammamix class that describes \deqn{g()  = pi0 delta1() + (1-pi0) gamma(shape, scale)} 
#' @param pi0 A scalar specifying the weight on point mass at 1
#' @param shape A scalar specifying the shape of gamma
#' @param scale A scalar specifying the scale of gamma
#'
#' @export 
point1_gamma <- function(pi0, shape, scale) {
  structure(data.frame(pi0 = pi0, shape = shape, scale = scale), class="point1_gamma")
}