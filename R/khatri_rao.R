#' khatri_rao computes the Khatri-Rao products between two matrices
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#' @return The Khatri-Rao product between x and y.
#' @title Khatri-Rao product
#' @noRd
khatri_rao <- function(x, y = NULL){
  if (is.null(x)) {
    return(y)
  }
  if (is.null(y)) {
    return(x)
  }
  r <- NCOL(x)
  vapply(
    seq_len(r),
    function(i) x[, i] %x% y[, i],
    FUN.VALUE = double(NROW(x) * NROW(y))
  )
}
