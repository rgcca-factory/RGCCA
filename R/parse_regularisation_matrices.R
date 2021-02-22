# Parse regularization matrices, invert them and take their square roots to make
# change of variables and solve the MGCCA optimization problem under easier
# constraints.
# @param reg_matrices If not NULL, a list. Each element
# of \eqn{reg_matrices} is a positive
# definite regularization matrix. There must be as many matrices as modes
# on the corresponding block and their dimensions must match the dimensions of
# the corresponding modes. 
# @param tau A regularization parameter. Value of 1 is always used for tensor
# blocks. If no numerical value is given for matrix block, a value is computed
# using the Schafer and Strimmer (2005) analytical formula.
# @param A An array that represents one block of variables Xj
# @param DIM A vector corresponding to the dimensions of array A.
# @param j An integer giving the number of the block A.
# @param bias A logical to know if a biased estimation of the covariance matrix
# should be considered.
# @return \item{P}{The same block of variables represented by \eqn{A} after 
# change of variables.}
# @return \item{M_inv_sqrt}{Either a matrix which represents the square 
# root of the inverse of the regularization matrix. If there are more than one
# regularization matrix, this is a list with each element being the list 
# described before for the different regularization matrices.}
# @return \item{tau}{The value of tau used to compute the regularization 
# matrix.}
# @title Parsing regularization matrices and applying change of variable
# @export parse_regularisation_matrices

parse_regularisation_matrices <- function(reg_matrices, tau, A, DIM, 
                                          j, bias = TRUE) {
  # TODO: add message if absolute values of eigenvalues are under a given tol
  sqrtMatrix = function(M, context = "matrix", d = NULL){
    eig        = eigen(M)
    if (any(abs(eig$values) < .Machine$double.eps)) {
      if (context == "matrix") {
        stop_rgcca(paste0("Regularized covariance matrix for block ", j, " is 
                          singular, try another value for tau[", j, "]."))
      } else {
        stop_rgcca(paste0("Mode ", d, " regularization matrix for block ", j, 
                          " is singular, please give an invertible matrix."))
      }
    }
    M_inv_sqrt = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    return(M_inv_sqrt)
  }
  P = (DIM[1]^(-1/2)) * A

  if (length(DIM) > 2){
    tau = 1
    if (is.null(reg_matrices)) {
      M_inv_sqrt = NULL
    } else if (is.list(reg_matrices)) {
      M_inv_sqrt = list()
      for (d in 1:length(reg_matrices)) {
        M_inv_sqrt[[d]] = sqrtMatrix(reg_matrices[[d]], context = "tensor", d=d)
        P               = mode_product(P, M_inv_sqrt[[d]], mode = d + 1)
      }
    }
    P = matrix(as.vector(P), nrow = DIM[1])
  }else{
    if (is.na(tau)) tau = tau.estimate(A)
    if(tau != 1){
      if (bias) {
        M = tau * diag(DIM[2]) + ((1 - tau)/(DIM[1] - 1)) * (t(A) %*% A)
      } else {
        M = tau * diag(DIM[2]) + ((1 - tau)/DIM[1]) * (t(A) %*% A)
      }
      M_inv_sqrt = sqrtMatrix(M, context = "matrix")
      P          = P %*% M_inv_sqrt
    }else{
      M_inv_sqrt = NULL
    }
  }
  return(list(P = P, M_inv_sqrt = M_inv_sqrt, tau = tau))
}
