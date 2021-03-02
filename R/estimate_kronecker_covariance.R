# estimate_kronecker_covariance estimates the covariance matrix of a set of
# random variables with an underlying tensor structure making the assumption
# that the real covariance matrix has a Kronecker structure.
# @param x A numerical array with at least 3 dimensions
# @return \item{factors}{List of estimated Kronecker factors.}
# @references Hoff, P. D. (2011), Separable covariance arrays via the Tucker
# product, with applications to multivariate relational data.
# Eun Jeong Min et al (2019), Tensor canonical correlation analysis.
# @title Kronecker covariance estimator
estimate_kronecker_covariance <- function(x) {
  DIM = dim(x)
  if (length(DIM) < 3) {
    stop_rgcca("x must have at least 3 dimensions")
  }
  factors = list()
  N       = DIM[1]
  D       = length(DIM) - 1
  x_bar   = apply(x, -1, mean)
  r       = (1 / N) * norm(matrix(x, nrow = N) - matrix(rep(x_bar, N), nrow = N, byrow = T), type = "F") ^ 2
  x_bar   = array(x_bar, dim = DIM[-1])
  for (d in 1:D) {
    x_bar_d      = unfold(x_bar, mode = d)
    x_bar_d      = x_bar_d %x% t(rep(1, N))
    x_d          = unfold(x, mode = d + 1)
    factors[[d]] = (1 / (N * r ^ ((D - 1) / D))) * tcrossprod(x_d - x_bar_d)
  }
  return(factors)
}

estimate_kronecker_mass <- function(x) {
  DIM = dim(x)
  N       = DIM[1]
  D       = length(DIM) - 1
  x_bar   = apply(x, -1, mean)
  r       = (1 / N) * norm(matrix(x, nrow = N) - matrix(rep(x_bar, N), nrow = N, byrow = T), type = "F") ^ 2
  return(r ^ ((D - 1) / D) )
}

generate_test_data <- function(DIM) {
  N       = DIM[1]
  factors = list()
  D       = length(DIM) - 1
  for (d in 1:D) {
    p                   = DIM[d + 1]
    Sigma               = diag(p)
    values              = runif(n = p * (p - 1) / 2)
    idx                 = matrix(1, p, p) - diag(p)
    idx[lower.tri(idx)] = 0
    Sigma[idx == 1]     = values / sqrt(p)
    Sigma               = Sigma + t(Sigma) - diag(p)
    factors[[d]]        = Sigma
  }
  Sigma = Reduce("%x%", rev(factors))
  x = mvrnorm(n = N, mu = rep(0, prod(DIM[-1])), Sigma = Sigma, empirical = F)
  return(list(x = array(x, dim = DIM), Sigma = Sigma))
}
