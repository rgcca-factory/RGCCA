#' Compute the R squared between a block and its extracted components.
#' @noRd
rsq <- function(y, X) {
  if (NCOL(X) == 1) {
    return(cor(X, y) ^ 2)
  }
  mean(apply(X, 2, function(x) {
    cor(x, y) ^ 2
  }))
}
