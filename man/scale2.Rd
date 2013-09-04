\name{scale2}
\alias{scale2}
\title{Scaling and Centering of Matrix-like Objects}
\usage{
  scale2(A, center = TRUE, scale = TRUE, bias = TRUE)
}
\arguments{
  \item{A}{A numeric matrix.}

  \item{center}{A logical value. If center = TRUE, each
  column is translated to have zero mean (default: TRUE).}

  \item{scale}{A logical value. If scale = TRUE, each
  column is transformed to have unit variance (default =
  TRUE).}

  \item{bias}{Logical value for biaised (\eqn{1/n}) or
  unbiaised (\eqn{1/(n-1)}) estimator of the var/cov
  (default = TRUE).}
}
\value{
  \item{A}{The centered and/or scaled matrix. The centering
  and scaling values (if any) are returned as attributes
  "scaled:center" and "scaled:scale".}
}
\description{
  Standardization (to zero means and unit variances) of
  matrix-like objects.
}

