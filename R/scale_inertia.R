#' Scale a list of blocks by dividing each block by its Frobenius norm
#' @inheritParams scale_lambda1
#' @noRd
scale_inertia <- function(blocks, sqrt_N, scale, na.rm) {
  blocks <- lapply(blocks, function(x) {
    if (na.rm) {
      z <- x
      z[is.na(z)] <- 0
    } else {
      z <- x
    }
    fac <- 1 / sqrt_N * norm(z, type = "F")
    y <- x / fac
    if (scale) {
      attr(y, "scaled:scale") <- attr(x, "scaled:scale") * fac
    } else {
      attr(y, "scaled:scale") <- rep(fac, NCOL(x))
    }
    return(y)
  })
}
