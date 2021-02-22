#' imputeMGCCA allows to choose the imputation method before running MGCCA   # TODO: ask if imputeMGCCA should not be mgccaNA
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams mgcca
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A symmetric matrix (J*J) that describes the relationships
#' between blocks}
#' @return \item{scheme}{A character or a function giving the link function for
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the squared values), "centroid" (the absolute values). Only, the horst
#'  scheme penalizes structural negative correlation. The factorial scheme
#'  discriminates more strongly the blocks than the centroid one}
#' @return \item{ncomp}{A vector of 1*J integers giving the number of component
#' for each blocks}
#' @return \item{crit}{A vector of integer that contains for each component the
#' values of the analysis criteria across iterations.}
#' @return \item{mode}{A \eqn{1 \times J} vector that contains the formulation
#' ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA
#' alogrithm}
#' @return \item{AVE}{A list of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized
#' Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A., Philippe C., and Frouin V. (2015). Kernel
#' Generalized Canonical Correlation Analysis. Computational Statistics and
#' Data Analysis, 90, 114-131.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statistical Applications in Genetics and Molecular Biology 4:32.
#' @title Multiway Generalized Canonical Correlation Analysis (MGCCA)
#' @export

mgccaNa=function(blocks, method, connection = 1 - diag(length(blocks)),
                 tau = rep(1, length(blocks)),
                 ncomp = rep(1, length(blocks)),
                 scheme = "centroid", scale = TRUE,
                 init = "svd", bias = TRUE, tol = 1e-08,
                 verbose = TRUE, scale_block = TRUE, prescaling = FALSE,
                 quiet = FALSE, regularisation_matrices = NULL, 
                 ranks = rep(1, length(blocks)))
{
	indNA=lapply(blocks, function(x){return(which(is.na(x), arr.ind = TRUE))})

  if(method=="complete") { A=intersection_list(blocks) }
  else if (is.function(method)) { A=method(blocks) }
	else { stop_rgcca("Only \"complete\" method is implemented to handle missing 
	                  data for MGCCA") }

  fit = mgcca(A, C = connection, tau = tau, ncomp = ncomp,
              verbose = verbose, scale = scale, init = init, bias = bias,
              scale_block = scale_block, scheme = scheme,
              tol = tol, prescaling = prescaling, quiet = quiet,
              regularisation_matrices = regularisation_matrices, 
              ranks = ranks)

 return(list(imputed_blocks = A, rgcca = fit, method, indNA = indNA))

}
