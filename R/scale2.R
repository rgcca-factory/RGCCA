#' Standardization (to zero means and unit variances) of matrix-like objects.
#' @inheritParams rgccad
#' @param A A numeric matrix.
#' @return The centered and scaled matrix. The centering and scaling
#' values (if any) are returned as attributes "scaled:center" and
#' "scaled:scale".
#' @title Scaling and Centering of Matrix-like Objects
#' @noRd
scale2 <- function(A, scale = TRUE, bias = TRUE) {
  # Center the data
  A <- scale(A, center = TRUE, scale = FALSE)

  # Scale if needed
  if (scale) {
    std <- sqrt(apply(A, 2, function(x) cov2(x, bias = bias)))
    std <- pmax(.Machine$double.eps, std) # Account for potentially 0 std
    A <- scale(A, center = FALSE, scale = std)
  }

  return(A)
}
