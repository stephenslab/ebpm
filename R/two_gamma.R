#' @title \code{two_gamma} class
#' @description two_gamma class that describes \deqn{g()  = pi0 gamma(shape1, scale1) + (1-pi0) gamma(shape2, scale2)} 
#' @param pi0 A scalar specifying the weight on point mass at 0
#' @param shape1,  A scalar specifying the shape of gamma in first component
#' @param scale1,  A scalar specifying the scale of gamma in first component
#' @param shape2,  A scalar specifying the shape of gamma in second component
#' @param scale2,  A scalar specifying the scale of gamma in second component
#'
#' @export 
two_gamma <- function(pi0, shape1, scale1, shape2, scale2) {
  structure(data.frame(pi0 = pi0, shape1 = shape1, scale1 = scale1, shape2 = shape2, scale2 = scale2), class="two_gamma")
}