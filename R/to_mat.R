#' Return a matrix version of a block. If it is already a matrix,
#' nothing is done, if it is an array, the first dimension is
#' kept while the others are concatenated and colnames are adapted.
#' @param block a block seen by RGCCA.
#' @return A matrix version of the block
#' @noRd
to_mat <- function(block) {
  if (is.matrix(block) || is.data.frame(block)) {
    return(block)
  }
  z <- matrix(block, nrow = NROW(block))
  rownames(z) <- rownames(block)
  grid <- expand.grid(dimnames(block)[-1])
  colnames(z) <- do.call(paste, c(grid, sep = "_"))
  return(z)
}
