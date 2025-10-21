#' Utility function to apply deflation
#' @param a List containing the computed weight vectors.
#' @param Y Matrix containing the computed components.
#' Each column corresponds to a block.
#' @param R List containing the deflated blocks.
#' @param P List containing the projection matrices used for deflation.
#' @param ndefl Vector of integers indicating the number of times each
#' block is deflated.
#' @param n Integer indicating the number of the current deflation.
#' @param superblock Logical indicating if there is a superblock.
#' @param comp_orth Logical indicating if the deflation leads to
#' orthogonal components.
#' @param response NULL or an integer indicating the position of a
#' response block.
#' @param na.rm Logical indicating if NA values should be removed.
#' @noRd
deflate <- function(a, Y, R, P, ndefl, n, superblock,
                    comp_orth, response, na.rm) {
  J <- length(a)
  pjs <- vapply(R, NCOL, FUN.VALUE = 1L)
  # Select the variable used to deflate blocks
  if (comp_orth) {
    var_defl <- lapply(seq_len(NCOL(Y)), function(i) Y[, i])
    left <- TRUE
  } else {
    var_defl <- a
    left <- FALSE
  }
  # If we aim for orthogonal components with a superblock, we need to deflate
  # the superblock and reconstruct the blocks from the superblock
  if (superblock && comp_orth) {
    defl_result <- deflation(R[[J]], var_defl[[J]], na.rm, left)
    P <- cbind(P, defl_result$p)
    cumsum_pjs <- cumsum(pjs)[seq_len(J - 1)]
    inf_pjs <- c(0, cumsum_pjs[seq_len(J - 2)]) + 1
    R <- lapply(seq(J - 1), function(b) {
      x <- defl_result$R[, seq(inf_pjs[b], cumsum_pjs[b]), drop = FALSE]
      colnames(x) <- colnames(defl_result$R)[seq(inf_pjs[b], cumsum_pjs[b])]
      return(x)
    })
    R[[J]] <- defl_result$R
  # Otherwise, the individual blocks are deflated
  } else {
    defl_result <- defl_select(var_defl, 
                               lapply(R, function(x) matrix(x, nrow = nrow(x))),
                               ndefl, n - 1, J,
                               na.rm = na.rm,
                               response = response,
                               left = left
    )
    R[seq_along(R)] <- lapply(seq_along(R), function(j) {
      x <- array(defl_result$resdefl[[j]], dim = dim(R[[j]]))
      dimnames(x) <- dimnames(R[[j]])
      return(x)
    })
    P <- lapply(seq(J), function(b) {
      cbind(P[[b]], defl_result$pdefl[[b]])
    })
    # If we aim for orthogonal weight vectors with a superblock, the superblock
    # must be reconstructed from the individual blocks
    if (superblock && !comp_orth) {
      R[[J]] <- do.call(cbind, R[seq(J - 1)])
    }
  }
  return(list(R = R, P = P))
}
