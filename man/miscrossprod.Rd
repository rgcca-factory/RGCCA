\name{miscrossprod}
\alias{miscrossprod}
\title{Cross product function for inputs with missing data.}
\usage{
  miscrossprod(x, y)
}
\arguments{
  \item{x}{A numeric vector.}

  \item{y}{A numeric vector.}
}
\value{
  \item{d.p}{The dot product between x and y: \eqn{x^ty}}
}
\description{
  Given vectors \eqn{x} and \eqn{y} as arguments, the
  function miscrossprod() returns the cross-product
  \eqn{x^ty}. miscrossprod() handles missing data.
}

