#' Scale a of block by dividing it by its first singular value
#' @inheritParams scaling
#' @param sqrt_N A numeric which corresponds to the square root of the number of
#' subjects, taking the bias into account.
#' @noRd
scale_lambda1 <- function(A, sqrt_N, scale, na.rm) {
  lambda <- sqrt(ifelse(
    NCOL(A) < NROW(A),
    eigen(
      pm(t(A / sqrt_N), A / sqrt_N, na.rm = na.rm),
      symmetric = TRUE, only.values = TRUE
    )$values[1],
    eigen(
      pm(A / sqrt_N, t(A / sqrt_N), na.rm = na.rm),
      symmetric = TRUE, only.values = TRUE
    )$values[1]
  ))
  if (lambda == 0) {
    return(A)
  }
  y <- A / lambda
  if (scale) {
    attr(y, "scaled:scale") <- attr(A, "scaled:scale") * lambda
  } else {
    attr(y, "scaled:scale") <- rep(lambda, NCOL(A))
  }
  return(y)
}
