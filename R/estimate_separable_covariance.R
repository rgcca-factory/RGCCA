#' estimate_separable_covariance estimates the covariance matrix of a set of
#' random variables with an underlying tensor structure making the assumption
#' that the real covariance matrix has a separable structure.
#' @param x A numerical array with at least 3 dimensions.
#' @return The list composed of the estimated terms in the separable covariance.
#' @references Hoff, P. D. (2011), Separable covariance arrays via the Tucker
#' product, with applications to multivariate relational data.
#' Eun Jeong Min et al (2019), Tensor canonical correlation analysis.
#' @title Separable covariance estimator
#' @noRd
estimate_separable_covariance <- function(x) {
  dim_x <- dim(x)
  n <- dim_x[1]
  d <- length(dim_x) - 1
  x_bar <- apply(x, -1, mean)
  r <- (1 / n) * norm(
    matrix(x, nrow = n) - matrix(rep(x_bar, n), nrow = n, byrow = TRUE),
    type = "F"
  )^2
  x_bar <- array(x_bar, dim = dim_x[-1])
  lapply(seq_len(d), function(m) {
    x_bar_m <- t(apply(x_bar, m, c))
    x_bar_m <- x_bar_m %x% t(rep(1, n))
    x_m <- t(apply(x, m + 1, c))
    (1 / (n * r^((d - 1) / d))) * tcrossprod(x_m - x_bar_m)
  })
}
