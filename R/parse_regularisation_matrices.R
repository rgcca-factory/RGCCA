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
# @return \item{P}{The same block of variables represented by \eqn{A} after 
# change of variables.}
# @return \item{M_inv}{Either a single list with values \eqn{M_sqrt} and 
# \eqn{Minv_sqrt} which represents respectively the square root and the square 
# root of the inverse of the regularization matrix. If there are more than one
# regularization matrix, this is a list with each element being the list 
# described before for the different regularization matrices.}
# @return \item{tau}{The value of tau used to compute the regularization 
# matrix.}
# @title Parsing regularization matrices and applying change of variable
# @export parse_regularisation_matrices

parse_regularisation_matrices <- function(reg_matrices, tau, A, DIM, bias = TRUE) {
  # TODO: add message if absolute values of eigenvalues are under a given tol
  sqrtMatrice = function(M){
    eig        = eigen(M)
    M_sqrt     = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    Minv_sqrt  = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    return(list(M_sqrt = M_sqrt, Minv_sqrt = Minv_sqrt))
  }
  P = (DIM[1]^(-1/2)) * A

  if (length(DIM) > 2){
    tau = 1
    if (is.null(reg_matrices)) {
      M_inv = NULL
    } else if (is.list(reg_matrices)) {
      M_inv = list()
      for (d in 1:length(reg_matrices)) {
        M_inv[[d]] = sqrtMatrice(reg_matrices[[d]])
        Minv_sqrt  = M_inv[[d]]$Minv_sqrt
        P          = mode_product(P, Minv_sqrt, mode = d + 1)
      }
    }
    P = matrix(as.vector(P), nrow = DIM[1])
  }else{
    if(tau != 1){
      if (!is.numeric(tau)) tau = tau.estimate(A)
      if (bias) {
        M = tau * diag(DIM[2]) + ((1 - tau)/(DIM[1]) - 1) * (t(A) %*% A)
      } else {
        M = tau * diag(DIM[2]) + ((1 - tau)/(DIM[1])) * (t(A) %*% A)
      }
      M_inv     = sqrtMatrice(M)
      Minv_sqrt = M_inv$Minv_sqrt
      P         = P %*% Minv_sqrt
    }else{
      M_inv = NULL
    }
  }
  return(list(P = P, M_inv = M_inv, tau = tau))
}
