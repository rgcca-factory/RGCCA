#' imputeRGCCA allows to choose the imputation method before running RGCCA
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams sgcca
#' @return \item{Y}{A list of J elements. Each element of Y is a
#' matrix that contains the RGCCA block components.}
#' @return \item{a}{A list of J elements. Each element of a is a
#' matrix that contains the outer weight vectors.}
#' @return \item{astar}{A list of J elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A symmetric matrix (J*J) that describes the relationships
#' between blocks.}
#' @return \item{scheme}{A character string defining the scheme function for
#' covariance maximization among "horst" (the identity function), "factorial"
#' (the squared values), "centroid" (the absolute values). The scheme argument
#' can be a user-defined function (e.g. function(x) x^4) but  have
#' to be a convex continuously differentiable function (to guaranty the
#' convergence properties of the algorithm).}
#' @return \item{ncomp}{A 1*J vector giving the number of component for each
#' block.}
#' @return \item{crit}{A numerical vector that contains for each component the
#' values of the RGCCA criteria across iterations.}
#' @return \item{mode}{A 1*J vector that contains the formulation
#' ("primal" or "dual") used for each block.}
#' @return \item{AVE}{A list of numerical values giving indicators of model
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
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA)
#' @export
#' @examples
#' data(Russett)
#' X_agric =as.matrix(Russett[, c("gini", "farm", "rent")])
#' X_ind = as.matrix(Russett[, c("gnpr", "labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' X_agric[c(2, 4), ] = NA # blockwise missing structure
#' X_ind[1, ] = NA # Ponctual NA
#' X_polit[5, 1] = NA
#' A = list(Agric = X_agric, Ind = X_ind, Polit = X_polit)
#' rgccaNa(A, method = "nipals")

sgccaNa=function(blocks, method, connection = 1 - diag(length(blocks)),
                 sparsity = rep(1, length(blocks)),
                 ncomp = rep(1, length(blocks)),
                 scheme = "centroid", scale = TRUE,
                 init = "svd", bias = TRUE, tol = 1e-08,
                 verbose = TRUE, scale_block = TRUE, prescaling = FALSE,
                 quiet = FALSE)
{
  shave.matlist <- function(mat_list, nb_cols)
    mapply(function(m,nbcomp) m[, 1:nbcomp, drop = FALSE],
           mat_list, nb_cols, SIMPLIFY = FALSE)

  shave.veclist <- function(vec_list, nb_elts)
	  mapply(function(m, nbcomp) m[1:nbcomp],
	         vec_list, nb_elts, SIMPLIFY = FALSE)
	indNA=lapply(blocks, function(x){return(which(is.na(x), arr.ind = TRUE))})
  na.rm=FALSE

  if(method=="complete") A=intersection_list(blocks)
  if(is.function(method)) A=method(blocks)
  if(method == "nipals") {na.rm = TRUE; A = blocks}

  fit = sgcca(A, C = connection, sparsity = sparsity, ncomp = ncomp,
              verbose = verbose, scale = scale,
              scale_block = scale_block, scheme = scheme,
              tol = tol, prescaling = prescaling, quiet = quiet)

 return(list(imputed_blocks = A, rgcca = fit, method, indNA = indNA))

}
