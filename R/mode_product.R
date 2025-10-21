#' Compute the mode product between a tensor X and a matrix y on a given mode m.
#' @param X An array.
#' @param y A matrix.
#' @param m A scalar designating a mode of the tensor X.
#' @examples
#' X <- array(rnorm(40 * 5 * 7), dim = c(40, 5, 7))
#' y <- rnorm(5)
#' res <- mode_product(X, y, m = 2)
#' print(dim(X))
#' print(dim(res))
#' @noRd
mode_product <- function(X, y, m = 1) {
  dim_X <- dim(X)

  # Unfold the tensor on dimension m
  perm <- c(m, seq_along(dim_X)[-m])
  X <- aperm(X, perm)

  # Compute the product
  X <- matrix(X, nrow = nrow(X))
  X <- t(y) %*% X

  # Reshape the result back to a tensor
  dim_X[m] <- NCOL(y)
  X <- array(X, dim = dim_X[perm])

  X <- aperm(X, c(1 + seq_len(m - 1), 1, m + seq_len(length(dim_X) - m)))
  return(X)
}
