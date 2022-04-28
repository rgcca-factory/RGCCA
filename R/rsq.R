#' Compute the R squared between a block and its extracted components.
#' @noRd
rsq <- function(y, X) {
  if (NCOL(X) == 1) {
    return(suppressWarnings(summary(lm(X ~ y))[["r.squared"]]))
  }
  mean(vapply(
    suppressWarnings(summary(lm(X ~ y))), "[[", "r.squared",
    FUN.VALUE = 1.0
  ))
}
