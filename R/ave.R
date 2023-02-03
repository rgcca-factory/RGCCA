#' Compute the AVE (Average Variance Explained) for each component for
#' a given block.
#'
#' @param X A matrix containing a given block.
#' @param Y A matrix containing the components associated to the given block.
#' @importFrom stats var
#' @noRd
ave <- function(X, Y) {
  var_X <- apply(X, 2, var, na.rm = TRUE)
  var_X <- var_X / sum(var_X)
  AVE_X <- drop(t(var_X) %*% cor(X, Y, use = "pairwise.complete.obs")^2)
  AVE_X_cum <- vapply(
    seq_len(NCOL(Y)),
    function(p) {
      Q <- qr.Q(qr(Y[, seq(p), drop = FALSE]))
      sum(t(var_X) %*% cor(X, Q, use = "pairwise.complete.obs")^2)
    },
    FUN.VALUE = 1.0
  )

  AVE_X_cor <- c(0, AVE_X_cum[-length(AVE_X_cum)])
  AVE_X_cor <- AVE_X_cum - AVE_X_cor
  return(list(AVE_X = AVE_X, AVE_X_cor = AVE_X_cor))
}
