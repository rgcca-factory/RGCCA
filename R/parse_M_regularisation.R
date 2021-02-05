# Parse regularization matrices, invert them and take their square roots to make
# change of variables and solve the MGCCA optimization problem under easier
# constraints.
# @param M_reg If not NULL, a list. Each element
# of \eqn{M_regn} is a symmetric positive
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
# @export parse_M_regularisation

parse_M_regularisation <- function(M_reg, tau, A, DIM) {

  sqrtMatrice = function(M){
    eig        = eigen(M)
    M_sqrt     = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    Minv_sqrt  = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    return(list(M_sqrt = M_sqrt, Minv_sqrt = Minv_sqrt))
  }
  P = (DIM[1]^(-1/2)) * A

  if (length(DIM) > 2){
    # We expect M_reg to be NULL or a list of matrices (one matrix per mode)
    tau = 1
    if (is.null(M_reg)) {
      M_inv = NULL
    } else if (is.list(M_reg)) {
      # Check dimensions and numbers of regularization matrices
      M_reg_DIM = lapply(M_reg, dim)
      if (any(sapply(M_reg_DIM, function(x) x[1]) != sapply(M_reg_DIM, function(x) x[2]))) {
        stop_rgcca("M_reg matrices must be square matrices")
      }
      if (length(M_reg) != length(DIM) - 1) {
        stop_rgcca("There should be as many M_reg matrices as modes")
      }
      if (any(sapply(M_reg_DIM, function(x) x[1]) != DIM[-1])) {
        stop_rgcca("M_reg matrices should match the mode dimensions")
      }
      M_inv = list()
      for (d in 1:length(M_reg)) {
        M_inv[[d]] = sqrtMatrice(M_reg[[d]])
        Minv_sqrt  = M_inv[[d]]$Minv_sqrt
        P          = mode_product(P, Minv_sqrt, mode = d + 1)
      }
    } else {
      stop_rgcca("For a tensor block, M_reg must be NULL or a list of matrices")
    }
    P = matrix(as.vector(P), nrow = DIM[1])
  }else{
    if(tau != 1){
      if (!is.numeric(tau)) tau = tau.estimate(A)
      M         = tau * diag(DIM[2]) + ((1 - tau)/(DIM[1])) * (t(A) %*% A)
      M_inv     = sqrtMatrice(M)
      Minv_sqrt = M_inv$Minv_sqrt
      P         = P %*% Minv_sqrt
    }else{
      M_inv = NULL
    }
  }
  return(list(P = P, M_inv = M_inv, tau = tau))
}
