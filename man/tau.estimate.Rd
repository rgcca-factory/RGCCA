\name{tau.estimate}
\alias{tau.estimate}
\title{Optimal shrinkage intensity parameters.}
\usage{
  tau.estimate(x)
}
\arguments{
  \item{x}{Data set on which the covariance matrix is
  estimated.}
}
\value{
  \item{tau}{Optimal shrinkage intensity parameter}
}
\description{
  Estimation of the optimal shrinkage parameters as
  described in [1,2] and implemented in
  \code{\link[SHIP:SHIP-package]{SHIP}} [2].
}
\references{
  [1] Schaefer J. and Strimmer K., 2005. A shrinkage
  approach to large-scale covariance matrix estimation and
  implications for functional genomics. Statist. Appl.
  Genet. Mol. Biol. 4:32.

  [2] Jelizarow M., Guillemot V., Tenenhaus A., Strimmer
  K., Boulesteix A.-L., 2010. Over-optimism in
  bioinformatics: an illustration. Bioinformatics
  26:1990-1998.
}

