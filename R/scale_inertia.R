#' Scale a block by dividing it by its Frobenius norm
#' @inheritParams scale_lambda1
#' @noRd
scale_inertia <- function(A, sqrt_N, scale, na.rm) {
  if (na.rm) {
    z <- A
    z[is.na(z)] <- 0
  } else {
    z <- A
  }
  fac <- 1 / sqrt_N * norm(z, type = "F")
  y <- A / fac
  if (scale) {
    attr(y, "scaled:scale") <- attr(A, "scaled:scale") * fac
  } else {
    attr(y, "scaled:scale") <- rep(fac, NCOL(A))
  }
  return(y)
}
