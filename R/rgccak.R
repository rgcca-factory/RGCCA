#' The function rgccak() is called by rgccad() and does not have to be used by
#' the user. The function rgccak() returns the RGCCA block components, block
#' weight vectors, etc., for each block and each dimension. Depending on the
#' dimensionality of each block \eqn{X_j, j = 1, ..., J}, the primal
#' (when n > p_j) or the dual (when n < p_j) algorithm is used (see
#' Tenenhaus et al. 2015)
#' @inheritParams rgccad
#' @param A  A list that contains the \eqn{J} blocks of variables from which
#' block components are constructed. It could be either the original matrices
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
#' @importFrom graphics plot
#' @importFrom Deriv Deriv
#' @noRd
rgccak <- function(A, C, tau = rep(1, length(A)), scheme = "centroid",
                   verbose = FALSE, init = "svd", bias = TRUE,
                   tol = 1e-08, na.rm = TRUE, n_iter_max = 1000) {
  if (is.function(scheme)) {
    g <- scheme
  } else {
    switch(scheme,
      "horst" = {
        g <- function(x) x
      },
      "factorial" = {
        g <- function(x) x^2
      },
      "centroid" = {
        g <- function(x) abs(x)
      }
    )
  }

  dg <- Deriv::Deriv(g, env = parent.frame())

  if (!is.numeric(tau)) {
    # From Schafer and Strimmer, 2005
    tau <- vapply(A, tau.estimate, na.rm = na.rm, FUN.VALUE = 1.0)
  }

  ### Initialization
  init_object <- rgcca_init(A, init, bias, na.rm, tau)
  a <- init_object$a
  Y <- init_object$Y

  iter <- 1
  crit <- NULL
  crit_old <- sum(C * g(cov2(Y, bias = bias)))
  a_old <- a

  repeat {
    update_object <- rgcca_update(A, bias, na.rm, tau, dg, C, a, Y, init_object)
    a <- update_object$a
    Y <- update_object$Y

    # Print out intermediate fit
    crit <- c(crit, sum(C * g(cov2(Y, bias = bias))))

    if (verbose) {
      cat(
        " Iter: ", formatC(iter, width = 3, format = "d"),
        " Fit: ", formatC(crit[iter], digits = 8, width = 10, format = "f"),
        " Dif: ", formatC(crit[iter] - crit_old,
          digits = 8, width = 10, format = "f"
        ), "\n"
      )
    }
    stopping_criteria <- c(
      drop(crossprod(unlist(a, FALSE, FALSE) - unlist(a_old, FALSE, FALSE))),
      abs(crit[iter] - crit_old)
    )

    if (any(stopping_criteria < tol) || (iter > 1000)) {
      break
    }

    crit_old <- crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  if (iter > n_iter_max) {
    warning(
      "The RGCCA algorithm did not converge after ", n_iter_max,
      " iterations."
    )
  }
  if (verbose) {
    if (iter <= n_iter_max) {
      message(
        "The RGCCA algorithm converged to a stationary point after ",
        iter - 1, " iterations \n"
      )
    }
    plot(crit, xlab = "iteration", ylab = "criteria")
  }

  result <- rgcca_postprocess(A, a, Y, g, na.rm)
  return(list(Y = result$Y, a = result$a, crit = crit, tau = tau))
}
