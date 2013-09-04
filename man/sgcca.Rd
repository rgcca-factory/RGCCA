\name{sgcca}
\alias{sgcca}
\title{Variable Selection For Generalized Canonical Correlation Analysis (SGCCA)}
\usage{
  sgcca(A, C = 1 - diag(length(A)), c1 = rep(1, length(A)),
    ncomp = rep(1, length(A)), scheme = "centroid",
    scale = TRUE, init = "svd", bias = TRUE,
    tol = .Machine$double.eps, verbose = FALSE)
}
\arguments{
  \item{A}{A list that contains the \eqn{J} blocks of
  variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ...,
  \mathbf{X_J}}.}

  \item{C}{A design matrix that describes the relationships
  between blocks (default: complete design).}

  \item{c1}{Either a \eqn{1 \times J} vector or a \eqn{\max
  (ncomp) \times J} matrix encoding the L1 constraints
  applied to the outer weight vectors. Elements of c1 vary
  between 0 and 1 (larger values of c1 correspond to less
  penalization). If c1 is a vector, L1-penalties are the
  same for all the weights corresponding to the same block
  but different components: \deqn{\forall h,
  \|a_{j,h}\|_{\ell_1} \leq c_1[j] \sqrt{p_j},} with
  \eqn{p_j} the number of variables of \eqn{\mathbf{X}_j}.
  If c1 is a matrix, each row \eqn{h} defines the
  constraints applied to the weights corresponding to
  components \eqn{h}: \deqn{\forall h, \|a_{j,h}\|_{\ell_1}
  \leq c_1[h,j] \sqrt{p_j}.}}

  \item{ncomp}{A \eqn{1 \times J} vector that contains the
  numbers of components for each block (default: rep(1,
  length(A)), which means one component per block).}

  \item{scheme}{Either "horst", "factorial" or "centroid"
  (Default: "centroid").}

  \item{scale}{If scale = TRUE, each block is standardized
  to zero means and unit variances (default: TRUE).}

  \item{init}{Mode of initialization use in the SGCCA
  algorithm, either by Singular Value Decompostion ("svd")
  or random ("random") (default : "svd").}

  \item{bias}{A logical value for biaised or unbiaised
  estimator of the var/cov.}

  \item{verbose}{Will report progress while computing if
  verbose = TRUE (default: TRUE).}

  \item{tol}{Stopping value for convergence.}
}
\value{
  \item{Y}{A list of \eqn{J} elements. Each element of Y is
  a matrix that contains the SGCCA components for each
  block.}

  \item{a}{A list of \eqn{J} elements. Each element of a is
  a matrix that contains the outer weight vectors for each
  block.}

  \item{astar}{A list of \eqn{J} elements. Each element of
  astar is a matrix defined as Y[[j]][, h] =
  A[[j]]\%*\%astar[[j]][, h]}

  \item{C}{A design matrix that describes the relationships
  between blocks (user specified).}

  \item{scheme}{The scheme chosen by the user (user
  specified).}

  \item{c1}{A vector or matrix that contains the value of
  c1 applied to each block \eqn{\mathbf{X}_j}, \eqn{ j=1,
  \ldots, J} and each dimension (user specified).}

  \item{ncomp}{A \eqn{1 \times J} vector that contains the
  number of components for each block (user specified).}

  \item{crit}{A vector that contains the values of the
  objective function at each iterations.}

  \item{AVE}{Indicators of model quality based on the
  Average Variance Explained (AVE): AVE(for one block),
  AVE(outer model), AVE(inner model).}
}
\description{
  SGCCA extends RGCCA to address the issue of variable
  selection. Specifically, RGCCA is combined with an
  L1-penalty that gives rise to Sparse GCCA (SGCCA) which
  is implemented in the function sgcca(). Given \eqn{J}
  matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ...,
  \mathbf{X_J}}, that represent \eqn{J} sets of variables
  observed on the same set of \eqn{n} individuals. The
  matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ...,
  \mathbf{X_J}} must have the same number of rows, but may
  (and usually will) have different numbers of columns.
  Bblocks are not necessarily fully connected within the
  SGCCA framework. Hence the use of SGCCA requires the
  construction (user specified) of a design matrix
  (\eqn{\mathbf{C}}) that characterizes the connections
  between blocks. Elements of the (symmetric) design matrix
  \eqn{\mathbf{C} = (c_{jk})} areequal to 1 if block
  \eqn{j} and block \eqn{k} are connected, and 0 otherwise.
  Hence, the use of SGCCA requires the construction (user
  specified) of a design matrix (\eqn{\mathbf{C}}) which
  characterizes the connections between blocks. The SGCCA
  algorithm is very similar to the RGCCA algorithm and
  keeps the same monotone convergence properties (i.e. the
  bounded criteria to be maximized increases at each step
  of the iterative procedure). Moreover, using a deflation
  strategy, sgcca() enables computation of several SGCCA
  block components (specified by ncomp) for each block.
  Block components for each block are guaranteed to be
  orthogonal when using this deflation strategy. The
  so-called symmetric deflation is considered in this
  implementation, i.e. each block is deflated with respect
  to its own component. Moreover, we stress that the
  numbers of components per block could differ from one
  block to another.
}
\examples{
#############
# Example 1 #
#############
\dontrun{
# Download the dataset's package at http://biodev.cea.fr/sgcca/.
# --> gliomaData_0.4.tar.gz

require(gliomaData)
data(ge_cgh_locIGR)

A <- ge_cgh_locIGR$multiblocks
Loc <- factor(ge_cgh_locIGR$y) ; levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
tau = c(1, 1, 0)

# rgcca algorithm using the dual formulation for X1 and X2
# and the dual formulation for X3
A[[3]] = A[[3]][, -3]
result.rgcca = rgcca(A, C, tau, ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE)
# sgcca algorithm
result.sgcca = sgcca(A, C, c1 = c(.071,.2, 1), ncomp = c(2, 2, 1),
                     scheme = "centroid", verbose = FALSE)

############################
# plot(y1, y2) for (RGCCA) #
############################
layout(t(1:2))
plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], col = "white", xlab = "Y1 (GE)",
     ylab = "Y2 (CGH)", main = "Factorial plan of RGCCA")
text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], Loc, col = as.numeric(Loc), cex = .6)
plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], col = "white", xlab = "Y1 (GE)",
     ylab = "Y2 (GE)", main = "Factorial plan of RGCCA")
text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], Loc, col = as.numeric(Loc), cex = .6)

############################
# plot(y1, y2) for (SGCCA) #
############################
layout(t(1:2))
plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1], col = "white", xlab = "Y1 (GE)",
     ylab = "Y2 (CGH)", main = "Factorial plan of SGCCA")
text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1], Loc, col = as.numeric(Loc), cex = .6)

plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2], col = "white", xlab = "Y1 (GE)",
     ylab = "Y2 (GE)", main = "Factorial plan of SGCCA")
text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2], Loc, col = as.numeric(Loc), cex = .6)

# sgcca algorithm with multiple components and different L1 penalties for each components
# (-> c1 is a matrix)
init = "random"
result.sgcca = sgcca(A, C, c1 = matrix(c(.071,.2, 1, 0.06, 0.15, 1), nrow = 2, byrow = TRUE),
                     ncomp = c(2, 2, 1), scheme = "factorial", scale = TRUE, bias = TRUE,
                     init = init, verbose = FALSE)
# number of non zero elements per dimension
apply(result.sgcca$a[[1]], 2, function(x) sum(x!=0))
     #(-> 145 non zero elements for a11 and 107 non zero elements for a12)
apply(result.sgcca$a[[2]], 2, function(x) sum(x!=0))
     #(-> 85 non zero elements for a21 and 52 non zero elements for a22)
init = "svd"
result.sgcca = sgcca(A, C, c1 = matrix(c(.071,.2, 1, 0.06, 0.15, 1), nrow = 2, byrow = TRUE),
                     ncomp = c(2, 2, 1), scheme = "factorial", scale = TRUE, bias = TRUE,
                     init = init, verbose = FALSE)}
}
\references{
  Tenenhaus et al. Variable Selection For Generalized
  Canonical Correlation Analysis. 2013. Submitted to
  Biostatistics.
}

