#' L1-, L2- or LG- proximal mapping of a vector of numerics
#'
#' USAGE !!!!
#'
#' @param vec, vector of numeric value
#' @param lambda, proximal mapping parameter
#' @param grp, vector describing the groups
#'
#' @return the L1-, L2- or LG- proximal mapping of a vector of numerics
#' @export
#'
#' @examples
#' x <- c(-0.1, 1, 0.5)
#' lamb <- 0.5
#' g <- c(1, 1, 2)
#' proxL1(x, lamb) # = (0, 0.5, 0)
#' proxL2(x, lamb) # ~= (-0.1, 0.6, 0.3)
#' proxLG(x, lamb, c(1, 2, 3)) # = (0, 0.5, 0)
#' proxLG(x, lamb, c(1, 1, 1)) # ~= (-0.1, 0.6, 0.3)
#' proxLG(x, lamb, g) # ~= (-0.1, 0.5, 0)
NULL
#' @rdname prox
#' @export
proxL1 <- function(vec, lambda) {
  return(sign(vec)*pmax(0, abs(vec) - lambda))
}
#' @rdname prox
#' @export
proxL2 <- function(vec, lambda) {
  norm2x <- normL2(vec)
  if (norm2x < .Machine$double.eps) return(0*vec)
  return(vec * max(0, (norm2x - lambda) / norm2x))
}
#' @rdname prox
#' @export
proxLG <- function(vec, lambda, grp) {
  return(ave(vec, grp, FUN = function(xsub) proxL2(xsub, lambda)))
}
#' @rdname prox
#' @export
proxLG2 <- function(vec, lambda, grp) {
  zenormg <- ave(vec, grp, FUN = normL2)
  boolu <- (zenormg >= lambda) + 0
  shrinku <- (1 - lambda/zenormg) * vec
  return(boolu * shrinku)
}

