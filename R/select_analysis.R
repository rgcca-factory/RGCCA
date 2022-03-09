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
#' @param blocks List of blocks.
#' @param response Numerical value giving the position of the response block.
#' When the response argument is filled the supervised mode is automatically
#' activated.
#' @param connection Symmetric matrix (J*J) that describes the relationships
#' between blocks. Elements of the connection matrix must be positive; but
#' usually equal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0
#' otherwise.
#' @param penalty Vector of length J (or character string for 'optimal'
#' setting) indicating the values of the tuning parameters.
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

select_analysis <- function(blocks,
                            connection = 1 - diag(length(blocks)),
                            penalty = rep(1, length(blocks)),
                            ncomp = rep(1, length(blocks)),
                            scheme = "centroid",
                            superblock = TRUE,
                            method = "rgcca",
                            quiet = FALSE,
                            response = NULL) {
  call <- list(
    ncomp = ncomp, scheme = scheme, penalty = penalty,
    superblock = superblock, connection = connection
  )
  J <- length(blocks)
  c_all <- function(J) {
    return(matrix(1, J, J))
  }
  c_response <- function(J) {
    x <- matrix(0, J, J)
    x[, J] <- x[J, ] <- 1
    x[J, J] <- 0
    return(x)
  }
  c_pair <- function(J) {
    return(1 - diag(J))
  }

  switch(method,
    "rgcca" = {
    },
    "sgcca" = {
    },
    "pca" = {
      check_nblocks(blocks, "pca")
      ncomp <- max(ncomp)
      scheme <- "horst"
      penalty <- c(1, 1)
      superblock <- TRUE
      connection <- c_pair(2)
    },
    "spca" = {
      check_nblocks(blocks, "spca")
      ncomp <- max(ncomp)
      scheme <- "horst"
      penalty <- c(penalty[1], 1)
      superblock <- TRUE
      connection <- c_pair(2)
    },
    "pls" = {
      check_nblocks(blocks, "pls")
      scheme <- "horst"
      penalty <- c(1, 1)
      superblock <- FALSE
      connection <- c_pair(2)
    },
    "spls" = {
      check_nblocks(blocks, "spls")
      scheme <- "horst"
      superblock <- FALSE
      connection <- c_pair(2)
    },
    "cca" = {
      check_nblocks(blocks, "cca")
      scheme <- "horst"
      penalty <- c(0, 0)
      superblock <- FALSE
      connection <- c_pair(2)
    },
    "ifa" = {
      check_nblocks(blocks, "ifa")
      scheme <- "horst"
      penalty <- c(1, 1)
      superblock <- FALSE
      connection <- c_pair(2)
    },
    "ra" = {
      check_nblocks(blocks, "ra")
      scheme <- "horst"
      penalty <- c(1, 0)
      superblock <- FALSE
      connection <- c_pair(2)
    },
    "gcca" = {
      ncomp <- max(ncomp)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "maxvar" = {
      ncomp <- max(ncomp)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "maxvar-b" = {
      ncomp <- max(ncomp)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "maxvar-a" = {
      ncomp <- max(ncomp)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "mcoa" = {
      ncomp <- max(ncomp)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "cpca-1" = {
      ncomp <- max(ncomp)
      scheme <- "horst"
      penalty <- c(rep(1, J), 0)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "cpca-2" = {
      ncomp <- max(ncomp)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "cpca-4" = {
      ncomp <- max(ncomp)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "hpca" = {
      ncomp <- max(ncomp)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      superblock <- TRUE
      connection <- c_response(J + 1)
    },
    "maxbet-b" = {
      scheme <- "factorial"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "maxbet" = {
      scheme <- "horst"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "maxdiff-b" = {
      scheme <- "factorial"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_pair(J)
    },
    "maxdiff" = {
      scheme <- "horst"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_pair(J)
    },
    "sabscor" = {
      scheme <- "centroid"
      penalty <- rep(0, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "ssqcor" = {
      scheme <- "factorial"
      penalty <- rep(0, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "ssqcov-1" = {
      scheme <- "factorial"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "ssqcov-2" = {
      scheme <- "factorial"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_pair(J)
    },
    "ssqcov" = {
      scheme <- "factorial"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_pair(J)
    },
    "sumcor" = {
      scheme <- "horst"
      penalty <- rep(0, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "sumcov-1" = {
      scheme <- "horst"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "sumcov-2" = {
      scheme <- "horst"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_pair(J)
    },
    "sumcov" = {
      scheme <- "horst"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_pair(J)
    },
    "sabscov-1" = {
      scheme <- "centroid"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_all(J)
    },
    "sabscov-2" = {
      scheme <- "centroid"
      penalty <- rep(1, J)
      superblock <- FALSE
      connection <- c_pair(J)
    }
  )

  # Generate warnings if some parameters have been modified
  if (!quiet) {
    modified_parameters <- vapply(names(call), function(n) {
      if (!identical(call[[n]], get(n))) {
        return(n)
      }
      return("0")
    }, "0")
    modified_parameters <- modified_parameters[modified_parameters != "0"]
    n <- length(modified_parameters)
    if (n == 1) {
      warning(
        "Choice of method '", method, "' overwrote parameter '",
        modified_parameters, "'."
      )
    } else if (n > 1) {
      warning(
        "Choice of method '", method, "' overwrote parameters '",
        paste(modified_parameters, collapse = "', '"), "'."
      )
    }
  }

  if (method %in% c("rgcca", "sgcca") && superblock) {
    ncomp <- max(ncomp)
    penalty <- c(penalty, 1)[seq(J + 1)]
    connection <- c_response(J + 1)
  }

  return(list(
    scheme = scheme,
    penalty = penalty,
    ncomp = ncomp,
    connection = connection,
    superblock = superblock
  ))
}
