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

tau.estimate <- function(x, na.rm = TRUE) {
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE) {
    stop_rgcca("The data matrix must be numeric!")
  }
  p <- NCOL(x)
  if (p == 1) {
    return(1)
  }
  n <- NROW(x)
  corm <- cor(x, use = "pairwise.complete.obs")

  xs <- scale(x, center = TRUE, scale = TRUE)
  v <- (n / ((n - 1)^3)) * (pm(t(xs^2), xs^2, na.rm = na.rm) - 1 / n * (pm(t(xs), xs, na.rm = na.rm))^2)
  diag(v) <- 0
  I <- diag(NCOL(x))
  d <- (corm - I)^2

  tau <- (sum(v)) / sum(d)
  tau <- max(min(tau, 1), 0)
  return(tau)
}
