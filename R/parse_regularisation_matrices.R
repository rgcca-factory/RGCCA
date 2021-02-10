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

# TODO: Decide what to do with the checks, should they be there too, should they
# call existing checks, should there be no checks?
parse_regularisation_matrices <- function(reg_matrices, tau, A, DIM, bias = TRUE) {

  sqrtMatrice = function(M){
    eig        = eigen(M)
    M_sqrt     = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    Minv_sqrt  = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    return(list(M_sqrt = M_sqrt, Minv_sqrt = Minv_sqrt))
  }
  P = (DIM[1]^(-1/2)) * A

  if (length(DIM) > 2){
    # We expect reg_matrices to be NULL or a list of matrices (one matrix per mode)
    tau = 1
    if (is.null(reg_matrices)) {
      M_inv = NULL
    } else if (is.list(reg_matrices)) {
      # Check dimensions and numbers of regularization matrices
      reg_DIM = lapply(reg_matrices, dim)
      if (any(sapply(reg_DIM, function(x) x[1]) != sapply(reg_DIM, function(x) x[2]))) {
        stop_rgcca("reg_matrices matrices must be square matrices")
      }
      if (length(reg_matrices) != length(DIM) - 1) {
        stop_rgcca("There should be as many reg_matrices matrices as modes")
      }
      if (any(sapply(reg_DIM, function(x) x[1]) != DIM[-1])) {
        stop_rgcca("reg_matrices matrices should match the mode dimensions")
      }
      M_inv = list()
      for (d in 1:length(reg_matrices)) {
        M_inv[[d]] = sqrtMatrice(reg_matrices[[d]])
        Minv_sqrt  = M_inv[[d]]$Minv_sqrt
        P          = mode_product(P, Minv_sqrt, mode = d + 1)
      }
    } else {
      stop_rgcca("For a tensor block, reg_matrices must be NULL or a list of matrices")
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
