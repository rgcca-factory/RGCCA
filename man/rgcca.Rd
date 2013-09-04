\name{rgcca}
\alias{rgcca}
\title{Regularized Generalized Canonical Correlation Analysis (RGCCA)}
\usage{
  rgcca(A, C = 1 - diag(length(A)),
    tau = rep(1, length(A)), ncomp = rep(1, length(A)),
    scheme = "centroid", scale = TRUE, init = "svd",
    bias = TRUE, tol = .Machine$double.eps, verbose = TRUE)
}
\arguments{
  \item{A}{A list that contains the \eqn{J} blocks of
  variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ...,
  \mathbf{X_J}}.}

  \item{C}{A design matrix that describes the relationships
  between blocks (default: complete design).}

  \item{tau}{tau is either a \eqn{1 \times J} vector or a
  \eqn{\mathrm{max}(ncomp) \times J} matrix, and contains
  the values of the shrinkage parameters (default: tau = 1,
  for each block and each dimension). If tau = "optimal"
  the shrinkage paramaters are estimated for each block and
  each dimension using the Schafer and Strimmer (2005)
  analytical formula . If tau is a \eqn{1\times J} numeric
  vector, tau[j] is identical across the dimensions of
  block \eqn{\mathbf{X}_j}. If tau is a matrix, tau[k, j]
  is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th
  residual matrix for block \eqn{j})}

  \item{ncomp}{A \eqn{1 \times J} vector that contains the
  numbers of components for each block (default: rep(1,
  length(A)), which gives one component per block.)}

  \item{scheme}{The value is "horst", "factorial" or
  "centroid" (default: "centroid").}

  \item{scale}{If scale = TRUE, each block is standardized
  to zero means and unit variances (default: TRUE).}

  \item{verbose}{If verbose = TRUE, the progress will be
  report while computing (default: TRUE).}

  \item{init}{The mode of initialization to use in RGCCA
  algorithm. The alternatives are either by Singular Value
  Decompostion ("svd") or random ("random") (Default:
  "svd").}

  \item{bias}{A logical value for biaised or unbiaised
  estimator of the var/cov (default: bias = TRUE).}

  \item{tol}{The stopping value for convergence.}
}
\value{
  \item{Y}{A list of \eqn{J} elements. Each element of
  \eqn{Y} is a matrix that contains the RGCCA components
  for the corresponding block.}

  \item{a}{A list of \eqn{J} elements. Each element of
  \eqn{a} is a matrix that contains the outer weight
  vectors for each block.}

  \item{astar}{A list of \eqn{J} elements. Each element of
  astar is a matrix defined as Y[[j]][, h] =
  A[[j]]\%*\%astar[[j]][, h].}

  \item{C}{A design matrix that describes the relation
  between blocks (user specified).}

  \item{tau}{A vector or matrix that contains the values of
  the shrinkage parameters applied to each block and each
  dimension (user specified).}

  \item{scheme}{The scheme chosen by the user (user
  specified).}

  \item{ncomp}{A \eqn{1 \times J} vector that contains the
  numbers of components for each block (user specified).}

  \item{crit}{A vector that contains the values of the
  criteria across iterations.}

  \item{mode}{A \eqn{1 \times J} vector that contains the
  formulation ("primal" or "dual") applied to each of the
  \eqn{J} blocks within the RGCCA alogrithm}

  \item{AVE}{indicators of model quality based on the
  Average Variance Explained (AVE): AVE(for one block),
  AVE(outer model), AVE(inner model).}
}
\description{
  Regularized Generalized Canonical Correlation Analysis
  (RGCCA) is a generalization of regularized canonical
  correlation analysis to three or more sets of variables.
  Given \eqn{J} matrices \eqn{\mathbf{X_1}, \mathbf{X_2},
  ..., \mathbf{X_J}} that represent \eqn{J} sets of
  variables observed on the same set of \eqn{n}
  individuals. The matrices \eqn{\mathbf{X_1},
  \mathbf{X_2}, ..., \mathbf{X_J}} must have the same
  number of rows, but may (and usually will) have different
  numbers of columns. The aim of RGCCA is to study the
  relationships between these \eqn{J} blocks of variables.
  It constitutes a general framework for many multi-block
  data analysis methods. It combines the power of
  multi-block data analysis methods (maximization of well
  identified criteria) and the flexibility of PLS path
  modeling (the researcher decides which blocks are
  connected and which are not). Hence, the use of RGCCA
  requires the construction (user specified) of a design
  matrix, (\eqn{\mathbf{C}}), that characterize the
  connections between blocks. Elements of the (symmetric)
  design matrix \eqn{\mathbf{C} = (c_{jk})} is equal to 1
  if block \eqn{j} and block \eqn{k} are connected, and 0
  otherwise. The objective is to find a fixed point of the
  stationary equations related to the RGCCA optimization
  problem. The function rgcca() implements a monotonically
  convergent algorithm (i.e. the bounded criteria to be
  maximized increases at each step of the iterative
  procedure) that is very similar to the PLS algorithm
  proposed by Herman Wold. Moreover, depending on the
  dimensionality of each block \eqn{\mathbf{X}_j}, \eqn{j =
  1, \ldots, J}, the primal (when \eqn{n > p_j}) algorithm
  or the dual (when \eqn{n < p_j}) algorithm is used (see
  Tenenhaus et al. 2013). Moreover, by deflation strategy,
  rgcca() allow to compute several RGCCA block components
  (specified by ncomp) for each block. Block components of
  each block are guaranteed to be orthogonal with the use
  of the deflation. The so-called symmetric deflation is
  considered in this implementation, i.e. each block is
  deflated with respect to its own component. It should be
  noted that the numbers of components per block can differ
  from one block to another.
}
\examples{
#############
# Example 1 #
#############
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
A = list(X_agric, X_ind, X_polit)
#Define the design matrix (output = C)
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
result.rgcca = rgcca(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE)
lab = as.vector(apply(Russett[, 10:12], 1, which.max))
plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white",
     xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)")
text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], Russett[, 1], col = lab, cex = .7)

############################################
# Example 2: RGCCA and mutliple components #
############################################
############################
# plot(y1, y2) for (RGCCA) #
############################
result.rgcca = rgcca(A, C, tau = rep(1, 3), ncomp = c(2, 2, 1),
                     scheme = "factorial", verbose = TRUE)
plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], col = "white", xlab = "Y1 (GE)",
    ylab = "Y2 (CGH)", main = "Factorial plan of RGCCA")
text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], Russett[, 1], col = lab, cex = .6)
plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], col = "white", xlab = "Y1 (GE)",
    ylab = "Y2 (GE)", main = "Factorial plan of RGCCA")
text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], Russett[, 1], col = lab, cex = .6)

######################################
# example 3: RGCCA and leave one out #
######################################
Ytest = matrix(0, 47, 3)
result.rgcca = rgcca(A, C, tau = rep(1, 3), ncomp = rep(1, 3),
                     scheme = "factorial", verbose = TRUE)

for (i in 1:nrow(Russett)){
 B = lapply(A, function(x) x[-i, ])
 B = lapply(B, scale2)
 resB = rgcca(B, C, tau = rep(1, 3), scheme = "factorial", scale = FALSE, verbose = FALSE)
 #  look for potential conflicting sign among components within the loo loop.
 for (k in 1:length(B)){
   if (cor(result.rgcca$a[[k]], resB$a[[k]]) >= 0)
     resB$a[[k]] = resB$a[[k]] else resB$a[[k]] = -resB$a[[k]]
 }
 Btest =lapply(A, function(x) x[i, ])
 Btest[[1]]=(Btest[[1]]-attr(B[[1]],"scaled:center"))/(attr(B[[1]],"scaled:scale"))
 Btest[[2]]=(Btest[[2]]-attr(B[[2]],"scaled:center"))/(attr(B[[2]],"scaled:scale"))
 Btest[[3]]=(Btest[[3]]-attr(B[[3]],"scaled:center"))/(attr(B[[3]],"scaled:scale"))
 Ytest[i, 1] = Btest[[1]]\%*\%resB$a[[1]]
 Ytest[i, 2] = Btest[[2]]\%*\%resB$a[[2]]
 Ytest[i, 3] = Btest[[3]]\%*\%resB$a[[3]]
}
lab = apply(Russett[, 10:12], 1, which.max)
plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white",
     xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Ind. Development)")
text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], Russett[, 1], col = lab)
text(Ytest[, 1], Ytest[, 2], substr(Russett[, 1], 1, 1), col = lab)
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

