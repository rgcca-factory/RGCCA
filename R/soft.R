#' the function soft() encodes the soft-thresholding operator
#' @param x numeric vector on which soft-thresholding is applied
#' @param d numeric, threshold
#' @noRd
soft <- function(x, d) {
  return(sign(x) * pmax(0, abs(x) - d))
}
