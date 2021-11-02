# Standardization (to zero means and unit variances) of matrix-like objects.
# @inheritParams rgccad
# @param A A numeric matrix.
# @param center A logical value. If center = TRUE, each column has zero mean.
# @return \item{A}{The centered and/or scaled matrix. The centering and scaling values (if any) are returned as attributes "scaled:center" and "scaled:scale".}
# @title Scaling and Centering of Matrix-like Objects
# @export scale2

scale2<-function (A, center = TRUE, scale = TRUE, bias = TRUE)
{
  if (center == TRUE & scale == TRUE) {
    A = scale(A, center = TRUE, scale = FALSE)
    std = sqrt(apply(A, 2, function(x) cov2(x, bias = bias)))
    if (any(std == 0)) {
      sprintf("there were %d constant variables", sum(std == 0))
    }
    A = sweep(A, 2, std, FUN = "/")
    attr(A, "scaled:scale") = std
    return(A)
  }

  if (center == TRUE & scale == FALSE) {
    A = scale(A, center = TRUE, scale = FALSE)
    return(A)
  }

  if (center == FALSE & scale == FALSE) {
    return(A)
  }

  if (center == FALSE & scale == TRUE) {
    std = apply(A, 2, function(x) cov2(x, bias = bias))
    A = sweep(A, 2, std, FUN = "/")
    attr(A, "scaled:scale") = std
    return(A)
  }
}
