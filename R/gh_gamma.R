#' @title \code{gh_gamma} class
#' @description gh_gamma class 
#' @param alpha 
#' @param a
#' @param b
#' @param phi
#' @param gam
#'
#' @export 
point_gamma <- function(alpha, a, b, phi, gam) {
  structure(data.frame(alpha = alpha, a = a, b = b, phi = phi, gam = gam), class="gh_gamma")
}