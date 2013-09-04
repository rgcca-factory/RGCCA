\name{sgccak}
\alias{sgccak}
\title{Internal function for computing the SGCCA parameters (SGCCA block components, outer weight vectors etc.)}
\usage{
  sgccak(A, C, c1 = rep(1, length(A)), scheme = "centroid",
    scale = FALSE, tol = .Machine$double.eps, init = "svd",
    bias = TRUE, verbose = TRUE)
}
\arguments{
  \item{A}{A list that contains the \eqn{J} blocks of
  variables from which block components are constructed. It
  could be eiher the original matrices (\eqn{\mathbf{X}_1,
  \mathbf{X}_2, \ldots, \mathbf{X}_J}) or the residual
  matrices (\eqn{\mathbf{X}_{h1}, \mathbf{X}_{h2}, \ldots,
  \mathbf{X}_{hJ}}).}

  \item{C}{A design matrix that describes the relationships
  between blocks.}

  \item{c1}{A \eqn{1 \times J} vector that contains the
  value of c1 applied to each block. The L1 bound on a[[j]]
  is \deqn{ \|a_{j}\|_{\ell_1} \leq c_1[j] \sqrt{p_j}.}
  with \eqn{p_j} the number of variables of
  \eqn{\mathbf{X}_j} and with c1[j] between 0 and 1 (larger
  L1 bound corresponds to less penalization).}

  \item{scheme}{Either "horst", "factorial" or "centroid"
  (default: centroid).}

  \item{scale}{If scale = TRUE, each block is standardized
  to zero means and unit variances (default: TRUE).}

  \item{init}{Mode of initialization of the SGCCA
  algorithm. Either by Singular Value Decompostion ("svd")
  or random ("random") (default: "svd").}

  \item{bias}{Logical value for biaised (\eqn{1/n}) or
  unbiaised (\eqn{1/(n-1)}) estimator of the var/cov.}

  \item{verbose}{Reports progress while computing, if
  verbose = TRUE (default: TRUE).}

  \item{tol}{Stopping value for convergence.}
}
\value{
  \item{Y}{A \eqn{n \times J} matrix of SGCCA block
  components.}

  \item{a}{A list of \eqn{J} elements. Each element
  contains the outer weight vector of each block.}

  \item{crit}{The values of the objective function at each
  iteration of the iterative procedure.}

  \item{converg}{Speed of convergence of the alogrithm to
  reach the tolerance.}

  \item{AVE}{Indicators of model quality based on the
  Average Variance Explained (AVE): AVE(for one block),
  AVE(outer model), AVE(inner model).}

  \item{C}{A design matrix that describes the relationships
  between blocks (user specified).}

  \item{scheme}{The scheme chosen by the user (user
  specified).}
}
\description{
  The function sgccak() is called by sgcca() and does not
  have to be used by the user. sgccak() enables the
  computation of SGCCA block components, outer weight
  vectors, etc., for each block and each dimension.
}

