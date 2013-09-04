\name{cov2}
\alias{cov2}
\title{Variance and Covariance (Matrices)}
\usage{
  cov2(x, y = NULL, bias = TRUE)
}
\arguments{
  \item{x}{A numeric vector, matrix or data.frame.}

  \item{y}{A numeric vector, matrix or data.frame.}

  \item{bias}{A logical value. If bias = TRUE, \eqn{n} is
  used to give a biased estimator of the (co)variance. If
  bias = FALSE, \eqn{n-1} is used (default: TRUE).}
}
\value{
  \item{C}{Estimation of the variance (resp. covariance) of
  x (resp. x and y).}
}
\description{
  cov2() is similar to cov() but has an additional
  argument. The denominator \eqn{n} (bias = TRUE) can be
  used (instead of \eqn{n-1}) to give a biased estimator of
  the (co)variance.
}

