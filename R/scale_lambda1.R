#' Scale a list of blocks by dividing each block by its first singular value
#' @inheritParams scale
#' @param sqrt_N A numeric which corresponds to the square root of the number of
#' subjects, taking the bias into account.
#' @noRd
scale_lambda1 <- function(blocks, sqrt_N, scale, na.rm) {
  blocks <- lapply(blocks, function(x) {
    lambda <- sqrt(ifelse(
      NCOL(x) < NROW(x),
      eigen(pm(t(x / sqrt_N), x / sqrt_N, na.rm = na.rm))$values[1],
      eigen(pm(x / sqrt_N, t(x / sqrt_N), na.rm = na.rm))$values[1]
    ))
    y <- x / lambda
    if (scale) {
      attr(y, "scaled:scale") <- attr(x, "scaled:scale") * lambda
    } else {
      attr(y, "scaled:scale") <- rep(lambda, NCOL(x))
    }
    return(y)
  })
}
