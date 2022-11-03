#' Compute the R squared between a block and its extracted components.
#' @noRd
rsq <- function(y, X) {
  if (NCOL(X) == 1) {
    return(cor(X, y) ^ 2)
  }
  reg <- suppressWarnings(lm(X ~ y))
  mean(
    1 - apply(reg$residuals, 2, function(x) sum(x^2)) / apply(X, 2, function(x) sum(x^2))
  )
}
