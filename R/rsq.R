#' Compute the R squared between a block and its extracted components.
#' @noRd
rsq <- function(y, X) {
  if (NCOL(y) == 1) {
    return(mean(cor(X, y) ^ 2))
  }
  reg <- suppressWarnings(lm(X ~ y))
  mean(
    1 - apply(
      matrix(reg$residuals, ncol = NCOL(reg$residuals)), 2, function(x) sum(x^2)
    ) / apply(X, 2, function(x) sum(x^2))
  )
}
