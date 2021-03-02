# Estimation of the optimal shrinkage parameters as described in [1,2] and
# implemented in a more general version within the SHIP package [2].
# @param x  Data set on which the covariance matrix is estimated.
# @param na.rm if TRUE, calculations are made on available data
# @return \item{tau}{Optimal shrinkage intensity parameter}
# @title Optimal shrinkage intensity parameters.
# @references [1] Schaefer J. and Strimmer K., 2005. A shrinkage approach to
# large-scale covariance matrix estimation and implications for functional
# genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
# @references [2] Jelizarow M., Guillemot V., Tenenhaus A., Strimmer K.,
# Boulesteix A.-L., 2010. Over-optimism in bioinformatics: an illustration.
# Bioinformatics 26:1990-1998.

tau.estimate <- function(x, na.rm=TRUE) {
  if (is.null(dim(x)) || (length(dim(x)) == 2 && dim(x)[2] == 1)) return(1)
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE)
    stop_rgcca("The data matrix must be numeric!")
  p <- NCOL(x)
  n <- NROW(x)
  corm <- cor(x,use="pairwise.complete.obs")
  if(na.rm)
  {
    nmat=t(!is.na(x))%*%(!is.na(x))
    nmat[nmat==0]=NA
    xs <- scale3(x, center = TRUE, scale = TRUE, bias = FALSE)
    v <- (nmat/((nmat - 1)^3)) * (pm(t(xs^2),xs^2) - 1/nmat * (pm(t(xs),xs))^2)
    diag(v) <- 0
    m <- matrix(rep(apply(pm(t(xs),xs), 2, mean), p), p, p)
    I <- diag(NCOL(x))
    d <- (corm - I)^2
  }
  else
  {
    xs <- scale(x, center = TRUE, scale = TRUE)
    v <- (n/((n - 1)^3)) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
    diag(v) <- 0
    m <- matrix(rep(apply(xs^2, 2, mean), p), p, p)
    I <- diag(NCOL(x))
    d <- (corm - I)^2
  }

  tau <- (sum(v))/sum(d)
  tau <- max(min(tau, 1), 0)
  return(tau)
}

tau.tensor.estimate <- function(x) {
  DIM   = dim(x)
  tau   = rep(1, length(DIM) - 1)
  N     = DIM[1]
  for (d in 1:(length(DIM) - 1)) {
    x_d     = unfold(x, mode = d + 1)
    S       = tcrossprod(x_d) / (N - 1)
    U       = diag(DIM[d + 1])
    tmp     = lapply(1:N, function(n) tcrossprod(unfold(x[n, ,], mode = d)) ^ 2)
    tmp     = Reduce("+", tmp)
    Var.S   = (tmp - tcrossprod(x_d) ^ 2 / N) * N * (N - 1) ^ (-3)
    diag(Var.S) = 0
    tau[d]      = sum(Var.S) / sum((S - U * S) ^ 2)
    # tau[d]      = sum(Var.S) / sum((S - U) ^ 2)
  }
  return(tau)
}
