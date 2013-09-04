\name{defl.select}
\alias{defl.select}
\title{deflation function}
\usage{
  defl.select(yy, rr, nncomp, nn, nbloc)
}
\arguments{
  \item{yy}{A matrix that contains the SGCCA block
  components of each block: \eqn{\mathbf{y}_{1h}, \ldots,
  \mathbf{y}_{Jh}}}

  \item{rr}{A list that contains the residual matrices
  \eqn{\mathbf{X}_{1h}, \ldots, \mathbf{X}_{Jh}}}

  \item{nncomp}{A \eqn{1 \times J} vector that contains the
  number of components to compute for each block.}

  \item{nn}{A \eqn{1 \times J} vector that contains the
  numbers of already computed components for each block}

  \item{nbloc}{Number of blocks.}
}
\value{
  \item{resdefl}{A list of \eqn{J} elements that contains
  \eqn{\mathbf{X}_{1,h+1}, \ldots, \mathbf{X}_{J,h+1}}.}

  \item{pdefl}{A list of \eqn{J} elements that contains
  \eqn{\mathbf{p}_{1h}, \ldots, \mathbf{p}_{Jh}}.}
}
\description{
  The function defl.select() computes residual matrices
  \eqn{\mathbf{X}_{1,h+1}, \ldots, \mathbf{X}_{J,h+1}}.
  These residual matrices are determined according to the
  following formula: \eqn{\mathbf{X}_{j,h+1} =
  \mathbf{X}_{jh} - \mathbf{y}_{jh} \mathbf{p}_{jh}^t}.
}

