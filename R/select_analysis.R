#' Define the parameters associated with each multi-block component
#' method of the literature.
#'
#' @param method A character string indicating the multi-block component
#' method to consider: rgcca, sgcca, pca, spca, pls, spls, cca,
#' ifa, ra, gcca, maxvar, maxvar-b, maxvar-a, mcoa,cpca-1, cpca-2,
#' cpca-4, hpca, maxbet-b, maxbet, maxdiff-b, maxdiff,
#' sabscor, ssqcor, ssqcov-1, ssqcov-2, ssqcov, sumcor,
#' sumcov-1, sumcov-2, sumcov, sabscov-1, sabscov-2.
#' @inheritParams plot_var_2D
#' @inheritParams set_connection
#' @param blocks A list that contains the J blocks of variables
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}{X1, X2, ..., XJ}.
#' Block \eqn{\mathbf{X}_j}{Xj} is a matrix of dimension
#' \eqn{n \times p_j}{n x p_j} where n is the number of
#' observations and \eqn{p_j} the number of variables.
#' @param response Numerical value giving the position of the response block.
#' When the response argument is filled the supervised mode is automatically
#' activated.
#' @param connection  A symmetric matrix (\eqn{J \times J}{J x J}) that
#' describes the relationships between blocks.
#' @param tau Either a \eqn{1 \times J}{1 x J} vector or a
#' \eqn{\mathrm{max}(ncomp) \times J}{max(ncomp) x J} matrix containing
#' the values of the regularization parameters (default: tau = 1, for each
#' block and each dimension). The regularization parameters varies from 0
#' (maximizing the correlation) to 1 (maximizing the covariance). If
#' tau = "optimal" the regularization parameters are estimated for each block
#' and each dimension using the Schafer and Strimmer (2005) analytical formula.
#' If tau is a \eqn{1 \times J}{1 x J} vector, tau[j] is identical across the
#' dimensions of block \eqn{\mathbf{X}_j}{Xj}. If tau is a matrix, tau[k, j]
#' is associated with \eqn{\mathbf{X}_{jk}}{Xjk} (kth residual matrix for
#' block j). The regularization parameters can also be estimated using
#' \link{rgcca_permutation} or \link{rgcca_cv}.
#' @param sparsity Either a \eqn{1*J} vector or a \eqn{max(ncomp) * J} matrix
#' encoding the L1 constraints applied to the outer weight vectors. The amount
#' of sparsity varies between \eqn{1/sqrt(p_j)} and 1 (larger values of sparsity
#' correspond to less penalization). If sparsity is a vector, L1-penalties are
#' the same for all the weights corresponding to the same block but different
#' components:
#' \deqn{for all h, |a_{j,h}|_{L_1} \le c_1[j] \sqrt{p_j},}
#' with \eqn{p_j} the number of variables of \eqn{X_j}.
#' If sparsity is a matrix, each row \eqn{h} defines the constraints applied to
#' the weights corresponding to components \eqn{h}:
#' \deqn{for all h, |a_{j,h}|_{L_1} \le c_1[h,j] \sqrt{p_j}.} It can be
#' estimated by using \link{rgcca_permutation}.
#' @param ncomp Vector of length J indicating the number of block components
#' for each block.
#' @param scheme Character string or a function giving the scheme function for
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the squared values), "centroid" (the absolute values). The scheme function
#'  can be any continously differentiable convex function and it is possible to
#'  design explicitely the scheme function (e.g. function(x) x^4) as argument of
#'  rgcca function.  See (Tenenhaus et al, 2017) for details.
#' @param quiet Logical value indicating if warning messages are reported.
#' @return \item{blocks}{List of blocks.}
#' @return \item{scheme}{Character string or a function giving the scheme
#' function used for covariance maximization.}
#' @return \item{penalty}{Vector of length J (or character string for
#' 'optimal' setting) indicating the values of the tuning parameters.}
#' @return \item{ncomp}{Vector of length J indicating the number of block
#' components or each block.}
#' @return \item{connection}{Symmetric matrix (J*J) that describes the
#' relationships between blocks.}
#' @return \item{superblock}{Logical value indicating if superblock is
#' included in the analysis.}
#' @return \item{response}{Integer giving the value of the response block if
#' any. NULL otherwise.}
#' @return \item{par}{String that indicates if penalty refers to tau or
#' sparsity.}
#' @return \item{gcca}{Function used to compute the analysis. Either rgccad
#' or sgcca.}

select_analysis <- function(blocks,
                            connection = 1 - diag(length(blocks)),
                            tau = rep(1, length(blocks)),
                            sparsity = rep(1, length(blocks)),
                            ncomp = rep(1, length(blocks)),
                            scheme = "centroid",
                            superblock = TRUE,
                            method = "rgcca",
                            quiet = FALSE,
                            response = NULL) {
  if (length(blocks) == 1) {
    if (sparsity == 1) {
      method <- "pca"
    } else {
      method <- "spca"
    }
  }

  method <- check_method(method)

  call <- list(
    ncomp = ncomp, scheme = scheme, tau = tau, sparsity = sparsity,
    superblock = superblock, connection = connection, response = response
  )
  J <- length(blocks)

  ### Utility functions to create the connection matrices
  name_c <- function(C, names_blocks) {
    rownames(C) <- colnames(C) <- names_blocks
    return(C)
  }
  c_all <- function(J, blocks) {
    return(name_c(matrix(1, J, J), names(blocks)))
  }
  c_response <- function(J, blocks, resp = J) {
    names_blocks <- names(blocks)
    if (J > length(blocks)) names_blocks <- c(names_blocks, "superblock")
    x <- matrix(0, J, J)
    x[, resp] <- x[resp, ] <- 1
    x[resp, resp] <- 0
    return(name_c(x, names_blocks))
  }
  c_pair <- function(J, blocks) {
    return(name_c(1 - diag(J), names(blocks)))
  }

  switch(method,
    "rgcca" = {
      par <- "tau"
      gcca <- rgccad
      penalty <- tau
    },
    "sgcca" = {
      par <- "sparsity"
      gcca <- sgcca
      penalty <- sparsity
    },
    "pca" = {
      check_nblocks(blocks, "pca")
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), 2)
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(2, blocks)
    },
    "spca" = {
      check_nblocks(blocks, "spca")
      par <- "sparsity"
      gcca <- sgcca
      ncomp <- rep(max(ncomp), 2)
      scheme <- "horst"
      penalty <- c(penalty[1], 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(2, blocks)
    },
    "pls" = {
      check_nblocks(blocks, "pls")
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- 2
      superblock <- FALSE
      connection <- c_response(2, blocks)
    },
    "spls" = {
      check_nblocks(blocks, "spls")
      par <- "sparsity"
      gcca <- sgcca
      scheme <- "horst"
      penalty <- check_penalty(sparsity, blocks, "sgcca")
      response <- 2
      superblock <- FALSE
      connection <- c_response(2, blocks)
    },
    "cca" = {
      check_nblocks(blocks, "cca")
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- c(0, 0)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(2, blocks)
    },
    "ifa" = {
      check_nblocks(blocks, "ifa")
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(2, blocks)
    },
    "ra" = {
      check_nblocks(blocks, "ra")
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- c(1, 0)
      response <- 2
      superblock <- FALSE
      connection <- c_response(2, blocks)
    },
    "gcca" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxvar" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxvar-b" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxvar-a" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "mcoa" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "cpca-1" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "horst"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "cpca-2" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "cpca-4" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "hpca" = {
      par <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxbet-b" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "maxbet" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "maxdiff-b" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "maxdiff" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sabscor" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "centroid"
      penalty <- rep(0, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "ssqcor" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(0, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "ssqcov-1" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "ssqcov-2" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "ssqcov" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sumcor" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(0, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sumcov-1" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "sumcov-2" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sumcov" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sabscov-1" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "centroid"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "sabscov-2" = {
      par <- "tau"
      gcca <- rgccad
      scheme <- "centroid"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    }
  )

  # Generate warnings if some parameters have been modified
  if (!quiet) {
    modified_parameters <- vapply(names(call), function(n) {
      if (n == par) assign(n, penalty)
      if (!identical(call[[n]], get(n))) {
        return(n)
      }
      return("0")
    }, character(1))
    modified_parameters <- modified_parameters[modified_parameters != "0"]
    warning(
      "Choice of method '", method, "' overwrote parameters '",
      paste(modified_parameters, collapse = "', '"), "'."
    )
  }

  if (method %in% c("rgcca", "sgcca")) {
    check_scheme(scheme)
    if (any(sparsity != 1)) {
      par <- "sparsity"
      gcca <- sgcca
      method <- "sgcca"
      penalty <- sparsity
    }
    if (!is.null(response)) {
      check_blockx("response", response, blocks)
      superblock <- FALSE
      connection <- c_response(J, blocks, resp = response)
    }
    if (superblock) {
      ncomp <- rep(max(ncomp), J + 1)
      connection <- c_response(J + 1, blocks)
      if (is.matrix(penalty)) {
        pen <- ifelse(ncol(penalty) < J + 1, 1, penalty[, J + 1])
        penalty <- cbind(penalty[, seq(J)], pen)
      } else {
        pen <- ifelse(length(penalty) < J + 1, 1, penalty[J + 1])
        penalty <- c(penalty[seq(J)], pen)
      }
    } else {
      connection <- check_connection(connection, blocks)
      ncomp <- check_ncomp(ncomp, blocks)
    }
    penalty <- check_penalty(penalty, blocks, method,
      superblock = superblock,
      ncomp = max(ncomp)
    )
  }

  return(list(
    scheme = scheme,
    penalty = penalty,
    ncomp = ncomp,
    connection = connection,
    superblock = superblock,
    response = response,
    gcca = gcca,
    par = par
  ))
}
