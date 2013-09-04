\name{rgccak}
\alias{rgccak}
\title{Internal function for computing the RGCCA parameters (RGCCA block components, outer weight vectors, etc.).}
\usage{
  rgccak(A, C, tau = "optimal", scheme = "centroid",
    scale = FALSE, verbose = FALSE, init = "svd",
    bias = TRUE, tol = .Machine$double.eps)
}
\arguments{
  \item{A}{A list that contains the \eqn{J} blocks of
  variables. Either the blocks (\eqn{\mathbf{X}_1,
  \mathbf{X}_2, \ldots, \mathbf{X}_J}) or the residual
  matrices (\eqn{\mathbf{X}_{h1}, \mathbf{X}_{h2}, \ldots,
  \mathbf{X}_{hJ}}).}

  \item{C}{A design matrix that describes the relationships
  between blocks. (Default: complete design).}

  \item{tau}{A \eqn{1 \times J} vector that contains the
  values of the shrinkage parameters \eqn{\tau_j}, \eqn{
  j=1, \ldots J}. (Default: \eqn{\tau_j = 1}, \eqn{ j=1,
  \ldots, J}). If tau = "optimal" the shrinkage intensity
  paramaters are estimated using the Schafer and Strimmer
  (2005) analytical formula.}

  \item{scheme}{Either "horst", "factorial" or "centroid"
  (default: centroid).}

  \item{scale}{if scale = TRUE, each block is standardized
  to zero means and unit variances (default: TRUE).}

  \item{verbose}{Will report progress while computing if
  verbose = TRUE (default: TRUE).}

  \item{init}{The mode of initialization to use in the
  RGCCA algorithm. The alternatives are either by Singular
  Value Decompostion or random (default : "svd").}

  \item{bias}{A logical value for either a biaised or
  unbiaised estimator of the var/cov.}

  \item{tol}{Stopping value for convergence.}
}
\value{
  \item{Y}{A \eqn{n \times J} matrix of RGCCA outer
  components}

  \item{Z}{A \eqn{n \times J} matrix of RGCCA inner
  components}

  \item{a}{A list of outer weight vectors}

  \item{crit}{The values of the objective function to be
  optimized in each iteration of the iterative procedure.}

  \item{converg}{Speed of convergence of the algorithm to
  reach the tolerance.}

  \item{AVE}{Indicators of model quality based on the
  Average Variance Explained (AVE): AVE(for one block),
  AVE(outer model), AVE(inner model).}

  \item{C}{A design matrix that describes the relationships
  between blocks (user specified).}

  \item{tau}{\eqn{1 \times J} vector containing the value
  for the tau penalties applied to each of the \eqn{J}
  blocks of data (user specified)}

  \item{scheme}{The scheme chosen by the user (user
  specified).}
}
\description{
  The function rgccak() is called by rgcca() and does not
  have to be used by the user. The function rgccak()
  computes the RGCCA block components, outer weight
  vectors, etc., for each block and each dimension.
  Depending on the dimensionality of each block
  \eqn{\mathbf{X}_j , j = 1, \ldots, J}, the primal (when
  \eqn{n > p_j}) or the dual (when \eqn{n < p_j}) algorithm
  is used (see Tenenhaus et al. 2013)
}
\references{
  Tenenhaus A. and Tenenhaus M., (2011), Regularized
  Generalized Canonical Correlation Analysis,
  Psychometrika, Vol. 76, Nr 2, pp 257-284.

  Tenenhaus A. et al., (2013), Kernel Generalized Canonical
  Correlation Analysis, submitted.

  Schafer J. and Strimmer K., (2005), A shrinkage approach
  to large-scale covariance matrix estimation and
  implications for functional genomics. Statist. Appl.
  Genet. Mol. Biol. 4:32.
}

