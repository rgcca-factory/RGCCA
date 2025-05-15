#' L1-norm, L2-norm or LG-norm of a vector of numerics
#'
#' USAGE !!!!
#'
#'
#' @param vec, vector of numeric value
#' @param grp, vector of numeric value
#'
#' @return the L1-norm, L2-norm or LG-norm of vec
#' @export
#'
#' @examples
#' x <- c(-0.1, 1, 0.5)
#' g <- c(1, 1, 2)
#' normL1(x) # = 1.6
#' normL2(x) # ~= 1.12
#' normLG(x, g) # ~= 1.5
#' normalizeL2(x)
NULL
#' @rdname norm
normL1 <- function(vec) {
  return(sum(abs(vec)))
}
#' @rdname norm
normL2 <- function(vec) {
  return(sqrt(sum(vec**2)))
}
#' @rdname norm
normLG <- function(vec, grp) {
  return(sum(tapply(vec, grp, normL2)))
}


