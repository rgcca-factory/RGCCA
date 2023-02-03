#' Compute the AVE (Average Variance Explained) for each component for
#' a given block.
#'
#' For PCA, the AVE of a component y associated to a block X is
#' ||Cov(X, y)||_2^2 / (Tr(Cov(X)) * Var(y)).
#' Since ||Cov(X, y)||_2^2 = Sum_{j = 1}^p Cov(X_j, y)^2, the AVE can be written
#' diag(Cov(X)) %*% Cor(X, y)^2 / Tr(Cov(X)).
#' The same formula is used for other methods. If components Y are not
#' orthogonal, a QR decomposition is sequentially applied and Y is replaced
#' by Q in the correlation.
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
