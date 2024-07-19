# workaround by Art Owen to avoid LAPACK errors
# See https://stat.ethz.ch/pipermail/r-help/2007-October/143508.html
svd_wrapper <- function(x, nu = min(n, p), nv = min(n, p), ...) {
  success <- FALSE
  n <- NROW(x)
  p <- NCOL(x)
  try({
    svd_x <- base::svd(x, nu, nv, ...)
    success <- TRUE
  }, silent = TRUE)
  if( success ) {
    return(svd_x)
  }
  try( {
    svd_tx <- base::svd(t(x), nv, nu, ...)
    success <- TRUE
  }, silent = TRUE )
  if( !success ) {
    stop("Error: svd(x) and svd(t(x)) both failed to converge.")
  }
  temp <- svd_tx$u
  svd_tx$u <- svd_tx$v
  svd_tx$v <- temp
  return(svd_tx)
}
