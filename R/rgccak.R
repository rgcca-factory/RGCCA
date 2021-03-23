#' The function rgccak() is called by rgccad() and does not have to be used by
#' the user. The function rgccak() returns the RGCCA block components, block
#' weight vectors, etc., for each block and each dimension. Depending on the
#' dimensionality of each block \eqn{X_j, j = 1, ..., J}, the primal
#' (when n > p_j) or the dual (when n < p_j) algorithm is used (see
#' Tenenhaus et al. 2015)
#' @inheritParams select_analysis
#' @inheritParams rgccad
#' @param A  A list that contains the \eqn{J} blocks of variables from which
#' block components are constructed. It could be eiher the original matrices
#' (\eqn{X_1, X_2, ..., X_J}) or the residual matrices
#' (\eqn{X_{h1}, X_{h2}, ..., X_{hJ}}).
#' @param C A symmetric matrix (J*J) that describes the
#' relationships between blocks.
#' @param na.rm If TRUE, RGCCA is run only on the available data (default value)
#' otherwise the NIPALS algorithm is used.
#' @return \item{Y}{Matrix of block components of dimension \eqn{n * J}}
#' @return \item{a}{A list of \eqn{J} elements. Each element of the list is a
#' matrix that contains a block weight vector associated with one block and
#' one deflation stage.}
#' @return \item{AVE_inner}{Average Variance Explained (AVE) of the inner model.}
#' @return \item{crit}{A list of max(ncomp) elements. Each element
#' (one per deflation stage) is a vector that contains the value of the RGCCA
#' objective function across iterations.}
#' @return \item{tau}{Either a 1*J vector or a \eqn{\mathrm{max}(ncomp)\times J}
#' matrix containing the values of the regularization parameters . The shrinkage
#' parameter tau varies from 0 (maximizing the correlation) to 1 (maximizing the
#' covariance). If tau = "optimal" the regularization paramaters are estimated
#' for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1\times J} vector, tau[j] is identical
#' across the dimensions of block \eqn{\mathbf{X}_j}. If tau is a matrix,
#' tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix
#' for block \eqn{j}). It can also be estimated by using
#' \link{rgcca_permutation}.}
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized
#' generalized canonical correlation analysis: A framework for sequential
#' multiblock component methods, Psychometrika, 82, 737-777
#' @references Tenenhaus A., Philippe C., and Frouin V. (2015). Kernel
#' Generalized Canonical Correlation Analysis. Computational Statistics and
#' Data Analysis, 90, 114-131.
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized
#' Canonical Correlation Analysis, Psychometrika, Vol. 76(2), pp 257-284.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statistical Applications in Genetics and Molecular Biology 4:32.
#' @title Internal function for computing the RGCCA parameters (RGCCA block
#' components, outer weight vectors, etc.).
#' @importFrom MASS ginv
#' @importFrom stats cor rnorm
#' @importFrom graphics plot
#' @importFrom Deriv Deriv

rgccak=function (A, C, tau = "optimal", scheme = "centroid", verbose = FALSE,
                 init = "svd", bias = TRUE, tol = 1e-08, na.rm = TRUE)
{

  if(mode(scheme) != "function")
  {
    if(!scheme %in% c("horst", "factorial", "centroid"))
      {stop_rgcca("Please choose scheme as 'horst', 'factorial', 'centroid'")}
    if(scheme == "horst"){ g <- function(x) x}
    if(scheme == "factorial"){ g <- function(x)  x^2}
    if(scheme == "centroid"){g <- function(x) abs(x)}
}
  else g <- scheme

    J <- length(A) # number of blocks
    n <- NROW(A[[1]]) # number of individuals
    pjs <- sapply(A, NCOL) # number of variables per block
    Y <- matrix(0, n, J)
    A <- lapply(A, as.matrix)
    a <- list()

    # Initialisation by SVD
    if (init == "svd") {
        for (j in 1:J) {
            a[[j]] <- initsvd(A[[j]], dual = FALSE)
        }
    }
    else if (init == "random") {
        for (j in 1:J) {
            a[[j]] <- rnorm(pjs[j]) # random initialisation
            a[[j]] <- a[[j]] / norm(a[[j]], type = "2")
        }
    }
    else {
        stop_rgcca("init should be either random or by SVD.")
    }

    for (j in 1:J) {
      Y[, j] <- pm(A[[j]] , a[[j]], na.rm = na.rm)
    }

    crit_old <- sum(C * g(cov2(Y, bias = bias)))
    iter = 1
    crit = numeric()
    Z = matrix(0, NROW(A[[1]]), J)
    a_old = a

    dg = Deriv::Deriv(g, env = parent.frame())

    repeat
    {
       for (j in 1:J)
      {
         dgx = dg(cov2(Y[, j], Y, bias = bias))
         Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) *
                            matrix(rep(dgx, n), n, J, byrow = TRUE) * Y)
         Az     = pm(t(A[[j]]), Z[, j], na.rm = TRUE)
	       a[[j]] = drop(1/sqrt(crossprod(Az))) * Az

		     Y[, j] = pm(A[[j]], a[[j]], na.rm = na.rm)
       }

      crit[iter] <- sum(C * g(cov2(Y, bias = bias)))
      if (verbose & (iter %% 1) == 0)
      {
          cat(" Iter: ", formatC(iter, width = 3, format = "d"),
              " Fit:", formatC(crit[iter], digits = 8,
                               width = 10, format = "f"),
              " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                                width = 10, format = "f"), "\n")
      }
       stopping_criteria = c(drop(crossprod(Reduce("c",
                                                   mapply("-", a, a_old)))),
                             crit[iter] - crit_old)

       if (any(stopping_criteria < tol) | (iter > 1000))
        {break}
      crit_old = crit[iter]
      a_old <- a
      iter <- iter + 1
    }
    if (iter > 1000)
        warning("The RGCCA algorithm did not converge after 1000 iterations.")
    if (iter < 1000 & verbose)
        cat("The RGCCA algorithm converged to a stationary point after",
            iter - 1, "iterations \n")
    if (verbose)
    {
        plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
    }
    AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)

    result <- list(Y = Y, a = a, crit = crit, AVE_inner = AVEinner, tau = tau)
    return(result)
}
