#' Compute the AVE (Average Variance Explained) for each component for
#' a given block.
#'
#' The AVE of a component is the variance of the component over the variance
#' of the associated block. In specific scenarios, components are not
#' orthogonal (comp_orth = FALSE, response block, superblock). Therefore,
#' they are sequentially orthogonalized to obtain cumulated AVEs. Individual
#' AVEs are then obtained by computing the increase of AVE due to the
#' additions of the components.
#' It can also happen that components are not in the span of the original
#' block (superblock), leading to AVE greater than 1. Such components are
#' projected on the span of the associated blocks before computing the AVE.
#' @param X A matrix containing a given block.
#' @param Y A matrix containing the components associated to the given block.
#' @param to_project A logical, if TRUE Y is projected to the span of X
#' before computing the AVE.
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
