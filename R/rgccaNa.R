#' imputeRGCCA allows to choose the imputation method before running RGCCA
#' @inheritParams select_analysis
#' @param blocks List that contains the J blocks of variables X1, X2, ..., XJ.
#' Block Xj is a matrix of dimension n x p_j where n is the number of
#' observations and p_j the number of variables.
#' @param method  Character string corresponding to the method used for
#' handling missing values ("nipals", "complete", "mean"). (default: "nipals").
#' \itemize{
#' \item{\code{"mean"}}{corresponds to imputation by the colmeans before
#' applying the S/RGCCA algorithm}
#' \item{\code{"complete"}}{corresponds to perform RGCCA on the fully observed
#' observations (observations with missing values are removed)}
#' \item{\code{"nipals"}}{corresponds to perform RGCCA algorithm on available
#' data (NIPALS-type algorithm)}}
#' @param tau Either a 1*J vector or a max(ncomp)*J matrix containing
#' the values of the regularization parameters (default: tau = 1, for each
#' block and each dimension). The regularization parameters varies from 0
#' (maximizing the correlation) to 1 (maximizing the covariance). If
#' tau = "optimal" the regularization paramaters are estimated for each block
#' and each dimension using the Schafer and Strimmer (2005) analytical formula .
#' If tau is a 1*J vector, tau[j] is identical across the dimensions
#' of block Xj. If tau is a matrix, tau[k, j] is associated with
#' X_jk (kth residual matrix for block j). The regularization parameters can
#' also be estimated using \link{rgcca_permutation} or \link{rgcca_cv}.
#' @param scale Logical value indicating if blocks are standardized.
#' @param scale_block Logical value indicating if each block is divided by
#' the square root of its number of variables.
#' @param verbose Logical value indicating if the progress of the
#' algorithm is reported while computing.
#' @param quiet Logical value indicating if warning messages are reported.
#' @param init Character string giving the type of initialization to use in
#' the  algorithm. It could be either by Singular Value Decompostion ("svd")
#' or by random initialisation ("random") (default: "svd").
#' @param bias A logical value for biaised (\eqn{1/n}) or unbiaised
#' (\eqn{1/(n-1)}) estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for the convergence of the algorithm.
#' @param prescaling Logical value indicating if the scaling has been done
#' outside of the function.
#' @return \item{Y}{List of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{List of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{List of \eqn{J} elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{Symmetric matrix (J*J) that describes the relationships
#' between blocks.}
#' @return \item{tau}{Regularization parameters used for the analysis.}
#' @return \item{ncomp}{Vector of 1*J integers giving the number of component
#' for each blocks.}
#' @return \item{crit}{Vector that contains for each component the
#' values of the objective function across iterations.}
#' @return \item{mode}{Vector of length J that contains the formulation
#' ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA
#' alogrithm.}
#' @return \item{AVE}{List of numeric values that reports indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized
#' Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A., Philippe C. and Frouin, V. (2015). Kernel
#' generalized canonical correlation analysis. Computational Statistics and
#' Data Analysis, 90, 114-131.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statistical Applications in Genetics and Molecular Biology 4:32.
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA)
#' @export
#' @examples
#' data(Russett)
#' X_agric = as.matrix(Russett[, c("gini", "farm", "rent")])
#' X_ind = as.matrix(Russett[,c("gnpr", "labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' X_agric[c(2,4),] = NA
#' X_ind[1, ]= NA
#' X_polit[5, 1]= NA
#' A = list(Agric = X_agric, Ind = X_ind, Polit = X_polit)
#' rgccaNa(A, method="nipals")
#'
#' blocks = list(agriculture = Russett[, seq(3)],
#'               industry = Russett[, 4:5],
#'               politic = Russett[, 6:11]
#'               )
#' blocks[[1]] = blocks[[1]][-c(3:4), ]
#' fit.rgcca = rgcca(blocks=blocks, type = "rgcca",
#'                   connection = 1-diag(3), scheme = "factorial",
#'                   tau = "optimal", ncomp = 2, verbose = TRUE,
#'                   superblock = TRUE)
#'
#'  # Define the label taking into account the blockwise missing structure
#'  lab = apply(fit.rgcca$A[[4]][, 9:11], 1, which.max)
#'  plot(fit.rgcca, type = "ind", resp = lab)
#'
rgccaNa=function (blocks, method, connection = 1 - diag(length(blocks)),
                  tau = rep(1, length(blocks)),
                  ncomp = rep(1, length(blocks)),
                  scheme = "centroid", scale = TRUE, init = "svd",
                  bias = TRUE, tol = 1e-08, verbose = TRUE,
                  scale_block=TRUE, prescaling = FALSE, quiet = FALSE)
{

    call=list(A = blocks, method = method, C = connection,
              tau = tau, ncomp = ncomp, scheme = scheme,
              scale = scale, init = init, bias = bias,
              tol = tol, verbose = verbose, scale_block = scale_block)

    shave.matlist <- function(mat_list, nb_cols)
      mapply(function(m,nbcomp) m[, 1:nbcomp, drop = FALSE],
             mat_list, nb_cols, SIMPLIFY = FALSE)
    shave.veclist <- function(vec_list, nb_elts)
      mapply(function(m, nbcomp) m[1:nbcomp],
             vec_list, nb_elts, SIMPLIFY = FALSE)

    indNA = lapply(blocks, function(x){return(which(is.na(x), arr.ind = TRUE))})
    na.rm = FALSE
    if(method == "complete") A = intersection_list(blocks)
    if(method == "nipals"){na.rm = TRUE ; A = blocks}

    fit = rgccad(A, C = connection, ncomp = ncomp, verbose = verbose,
                 scale = scale, init = init,
                 scale_block = scale_block,
                 tau = tau, scheme = scheme, tol = tol,
                 prescaling = prescaling, quiet = quiet)

    out = list(imputed_blocks = A, rgcca = fit, method, indNA = indNA)

    return(out)
}
