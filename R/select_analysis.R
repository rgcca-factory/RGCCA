#' Define the parameters associated with each multi-block component
#' method of the literature.
#'
#' @inheritParams rgcca
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
#' @return \item{param}{String that indicates if penalty refers to tau or
#' sparsity.}
#' @return \item{gcca}{Function used to compute the analysis. Either rgccad
#' or sgcca.}
#' @noRd
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
      param <- "tau"
      gcca <- rgccad
      penalty <- tau
    },
    "sgcca" = {
      param <- "sparsity"
      gcca <- sgcca
      penalty <- sparsity
    },
    "pca" = {
      check_nblocks(blocks, "pca")
      param <- "tau"
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
      param <- "sparsity"
      gcca <- sgcca
      ncomp <- rep(max(ncomp), 2)
      scheme <- "horst"
      penalty <- c(sparsity[1], sparsity[1])
      response <- NULL
      superblock <- TRUE
      connection <- c_response(2, blocks)
    },
    "pls" = {
      check_nblocks(blocks, "pls")
      param <- "tau"
      gcca <- rgccad
      ncomp[2] <- ncomp[1]
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- 2
      superblock <- FALSE
      connection <- c_response(2, blocks)
    },
    "spls" = {
      check_nblocks(blocks, "spls")
      param <- "sparsity"
      gcca <- sgcca
      ncomp[2] <- ncomp[1]
      scheme <- "horst"
      penalty <- check_penalty(sparsity, blocks, "sgcca")
      response <- 2
      superblock <- FALSE
      connection <- c_response(2, blocks)
    },
    "cca" = {
      check_nblocks(blocks, "cca")
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- c(0, 0)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(2, blocks)
    },
    "ifa" = {
      check_nblocks(blocks, "ifa")
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(2, blocks)
    },
    "ra" = {
      check_nblocks(blocks, "ra")
      param <- "tau"
      gcca <- rgccad
      ncomp[2] <- ncomp[1]
      scheme <- "horst"
      penalty <- c(1, 0)
      response <- 2
      superblock <- FALSE
      connection <- c_response(2, blocks)
    },
    "gcca" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxvar" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxvar-b" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxvar-a" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "mcoa" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "cpca-1" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "horst"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "cpca-2" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "cpca-4" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "hpca" = {
      param <- "tau"
      gcca <- rgccad
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      response <- NULL
      superblock <- TRUE
      connection <- c_response(J + 1, blocks)
    },
    "maxbet-b" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "maxbet" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "maxdiff-b" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "maxdiff" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sabscor" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "centroid"
      penalty <- rep(0, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "ssqcor" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(0, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "ssqcov-1" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "ssqcov-2" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "ssqcov" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sumcor" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(0, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sumcov-1" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "sumcov-2" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sumcov" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_pair(J, blocks)
    },
    "sabscov-1" = {
      param <- "tau"
      gcca <- rgccad
      scheme <- "centroid"
      penalty <- rep(1, J)
      response <- NULL
      superblock <- FALSE
      connection <- c_all(J, blocks)
    },
    "sabscov-2" = {
      param <- "tau"
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
      if (n == param) assign(n, penalty)
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
      param <- "sparsity"
      gcca <- sgcca
      method <- "sgcca"
      penalty <- sparsity
    }
    if (!is.null(response)) {
      check_blockx("response", response, blocks)
      ncomp[response] <- max(ncomp[-response])
      superblock <- FALSE
      connection <- c_response(J, blocks, resp = response)
    }
    if (superblock) {
      ncomp <- rep(max(ncomp), J + 1)
      connection <- c_response(J + 1, blocks)
      if (is.matrix(penalty)) {
        if (ncol(penalty) < J + 1) {
          pen <- 1
        } else {
          pen <- penalty[, J + 1]
        }
        penalty <- cbind(penalty[, seq(J)], pen)
      } else {
        pen <- ifelse(length(penalty) < J + 1, 1, penalty[J + 1])
        penalty <- c(penalty[seq(J)], pen)
      }
    } else {
      connection <- check_connection(connection, blocks)
    }
    penalty <- check_penalty(penalty, blocks, method,
      superblock = superblock,
      ncomp = max(ncomp)
    )
  }
  ncomp <- check_ncomp(
    ncomp, blocks,
    superblock = superblock, response = response
  )

  return(list(
    scheme = scheme,
    penalty = penalty,
    ncomp = ncomp,
    connection = connection,
    superblock = superblock,
    response = response,
    method = method,
    gcca = gcca,
    param = param
  ))
}
