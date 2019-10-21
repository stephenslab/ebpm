#' @param pi A vector of mixture proportions.
#' @param a A vector of shapes
#' @param b A vector of rates
#'
#' @export 
gammamix <- function(pi, a, b) {
structure(data.frame(pi, a, b), class="gammamix")
}

scale2gammamix_init <- function(scale){
  a = scale$a
  b = scale$b
  n = length(b)
  return(gammamix(pi = replicate(n, 1)/n, a = a, b = b))
}