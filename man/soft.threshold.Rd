\name{soft.threshold}
\alias{soft.threshold}
\title{The function soft.threshold() soft-thresholds a vector such that the L1-norm constraint is satisfied.}
\usage{
  soft.threshold(x, sumabs = 1)
}
\arguments{
  \item{x}{A numeric vector.}

  \item{sumabs}{A numeric constraint on x's L1 norm.}
}
\value{
  Returns a vector resulting from the soft thresholding of
  \eqn{x} given sumabs
}
\description{
  The function soft.threshold() soft-thresholds a vector
  such that the L1-norm constraint is satisfied.
}
\examples{
x <- rnorm(10)
soft.threshold(x,0.5)
}
\keyword{manip}

