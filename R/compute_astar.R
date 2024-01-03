#' Utility function to compute astar for a given GCCA model.
#' @param a List containing the computed weight vectors.
#' @param P List containing the projection matrices used for deflation.
#' @param superblock Logical indicating if there is a superblock.
#' @param comp_orth Logical indicating if the deflation leads to
#' orthogonal components.
#' @param N Integer indicating the number of times blocks are deflated.
#' @noRd
compute_astar <- function(a, P, superblock, supergroup, comp_orth, N) {
  J <- length(a)
  # If there is a superblock and components are orthogonal, astar is only
  # available for the superblock, same for the supergroup
  if ((superblock || supergroup) && comp_orth) {
    astar <- a[[J]]
    for (n in seq_len(N)) {
      astar[, n + 1] <- a[[J]][, n + 1] -
        astar[, seq(n), drop = FALSE] %*%
        drop(t(a[[J]][, n + 1]) %*% P[, seq(n), drop = FALSE])
    }
  } else {
    astar <- a
    # If weight vectors are orthogonal, astar is directly equal to a.
    if (comp_orth) {
      for (n in seq_len(N)) {
        astar <- lapply(seq(J), function(b) {
          cbind(
            astar[[b]][, seq(n), drop = FALSE],
            a[[b]][, n + 1] - astar[[b]][, seq(n), drop = FALSE] %*%
              drop(t(a[[b]][, n + 1]) %*% P[[b]][, seq(n), drop = FALSE])
          )
        })
      }
    }
  }
  return(astar)
}
