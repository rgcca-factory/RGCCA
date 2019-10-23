#'Uses cov3 to scale a list of matrix (taking into account the missing values)
#'
#'This method is used for the implementation of EM algorithm for missing data
#'
#' @param A A list of J blocks
#' @param center if TRUE, all variables are centered
#' @param scale if TRUE, all variables are scaled
#' @param bias if TRUE, the estimator of variance is SS/sqrt(n-1), if FALSE, it is SS/sqrt(n)
#' @return \item{A}{The resulting list of matrices}
#' @title scale3
#' @examples 
#'  data();...
scale3=function (A, center = TRUE, scale = TRUE, bias = TRUE) 
{
  if (center == TRUE & scale == TRUE) {
    A = scale(A, center = TRUE, scale = FALSE)
    std = sqrt(apply(A, 2, function(x) cov3(x, bias = bias)))
    if (any(std == 0)) {
      sprintf("there were %d constant variables", sum(std ==  0))
      std[std == 0] = 1
    }
    A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    return(A)
  }
  if (center == TRUE & scale == FALSE) {
    A = scale(A, center = TRUE, scale = FALSE)
    return(A)
  }
  if (center == FALSE & scale == TRUE) {
    std = apply(A, 2, function(x) cov3(x, bias = bias))
    A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    return(A)
  }
}