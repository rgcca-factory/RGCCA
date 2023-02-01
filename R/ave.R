#' Compute the AVE (Average Variance Explained) for each component for a given
#' block.
#' @noRd
ave <- function(X, Y, to_project) {
  # Replace missing values with 0 to be able to compute Frobenius norm
  if (any(is.na(X))) {
    X[is.na(X)] <- 0
  }

  # Project Y on X if to_project
  if (to_project) {
    Q <- qr.Q(qr(X))
    Y <- Q %*% (t(Q) %*% Y)
  }

  var_tot <- norm(X, type = "F")^2
  var_comp <- apply(Y, 2, norm, type = "2")^2
  AVE_X <- var_comp / var_tot
  AVE_X_cum <- vapply(
    seq_len(NCOL(Y)),
    function(p) {
      R <- qr.R(qr(Y[, seq(p), drop = FALSE]))
      sum(diag(R)^2) / var_tot
    },
    FUN.VALUE = 1.0
  )
  AVE_X_cor <- c(0, AVE_X_cum[-length(AVE_X_cum)])
  AVE_X_cor <- AVE_X_cum - AVE_X_cor
  return(list(AVE_X = AVE_X, AVE_X_cor = AVE_X_cor))
}
