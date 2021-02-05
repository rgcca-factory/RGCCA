# Computes the mode-k product between a tensor and a matrix or vector.
# @param A  A tensor.
# @param M A matrix or a vector.
# @param mode A mode to make the product along (default is 1).
# @param drop A boolean to know if dimension should be dropped after 
# multiplication by a vector.
# @return \item{A}{The result of the mode-k product}
# @title Tensor mode-k product

mode_product <- function(A, M, mode = 1, drop = FALSE) {
  # Exchange mode 1 with chosen mode
  idx = 1:length(dim(A))
  idx[mode] = 1
  idx[1] = mode
  A = aperm(A, idx)
  # Matricize the array and multiply by M
  DIM = dim(A)
  A = matrix(as.vector(A), nrow = dim(A)[1])
  M = matrix(M, nrow = NROW(M))
  A = t(M) %*% A
  DIM[1] = NCOL(M)
  # Fold back to array and exchange back modes
  A = array(A, dim = DIM)
  A = aperm(A, idx)
  if (drop) {
    A = drop(A)
  }
  return(A)
}
