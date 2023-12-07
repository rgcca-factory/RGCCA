#' Standardization (to zero means and unit variances) of matrix-like objects.
#' @inheritParams rgccad
#' @param A A numeric matrix.
#' @return The centered and scaled matrix. The centering and scaling
#' values (if any) are returned as attributes "scaled:center" and
#' "scaled:scale".
#' @title Scaling and Centering of Matrix-like Objects
#' @noRd
scale2_mg <- function(A, scale = TRUE, bias = TRUE, groups = NULL) {
  if (!is.null(groups)){
    A <- scale(A, center = TRUE, scale = FALSE)
    col_norm <- apply(A, 2, function(x) {norm(x, type = "2")})
    A <- scale(A, center = FALSE, scale = col_norm)
    return(A)
  }
  
  if (scale) {
    A <- scale(A, center = TRUE, scale = FALSE)
    std <- sqrt(apply(A, 2, function(x) cov2(x, bias = bias)))
    A <- scale(A, center = FALSE, scale = std)
    return(A)
  }
  
  A <- scale(A, center = TRUE, scale = FALSE)
  return(A)
}
