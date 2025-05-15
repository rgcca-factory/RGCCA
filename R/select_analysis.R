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
#' @importFrom utils modifyList
#' @noRd
select_analysis <- function(rgcca_args, blocks) {
  tau <- rgcca_args$tau
  ncomp <- rgcca_args$ncomp
  quiet <- rgcca_args$quiet
  scheme <- rgcca_args$scheme
  method <- rgcca_args$method
  response <- rgcca_args$response
  sparsity <- rgcca_args$sparsity
  comp_orth <- rgcca_args$comp_orth
  connection <- rgcca_args$connection
  superblock <- rgcca_args$superblock
  scale_block <- rgcca_args$scale_block
  lambda <- rgcca_args$lambda
  graph_laplacians <- rgcca_args$graph_laplacians

  if (length(blocks) == 1) {
    if (sparsity[1] == 1) {
      method <- "pca"
    } else {
      method <- "spca"
    }
  }

  method <- check_method(method)

  call <- list(
    ncomp = ncomp, scheme = scheme, tau = tau, sparsity = sparsity,
    superblock = superblock, connection = connection, response = response,
    comp_orth = comp_orth, scale_block = scale_block
  )
  J <- length(blocks)

  param_list <- list()

  switch(method,
    "rgcca" = {
      param <- "tau"
      penalty <- tau
    },
    "sgcca" = {
      param <- "sparsity"
      penalty <- sparsity
    },
    "netsgcca" = {
      param <- "sparsity"
      gcca <- netsgcca
      penalty <- sparsity
      lambda <- lambda
      graph_laplacians <- graph_laplacians
      param_list <- list(lambda = lambda, graph_laplacians = graph_laplacians)
    },
    "pca" = {
      check_nblocks(blocks, "pca")
      param <- "tau"
      ncomp <- rep(max(ncomp), 2)
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = 2)
    },
    "spca" = {
      check_nblocks(blocks, "spca")
      param <- "sparsity"
      ncomp <- rep(max(ncomp), 2)
      scheme <- "horst"
      penalty <- c(sparsity[1], sparsity[1])
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = 2)
    },
    "pls" = {
      check_nblocks(blocks, "pls")
      param <- "tau"
      ncomp[2] <- ncomp[1]
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- 2
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "response", J = 2)
    },
    "spls" = {
      check_nblocks(blocks, "spls")
      param <- "sparsity"
      ncomp[2] <- ncomp[1]
      scheme <- "horst"
      penalty <- check_penalty(sparsity, blocks, "sgcca")
      response <- 2
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "response", J = 2)
    },
    "cca" = {
      check_nblocks(blocks, "cca")
      param <- "tau"
      scheme <- "horst"
      penalty <- c(0, 0)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair", J = 2)
    },
    "ifa" = {
      check_nblocks(blocks, "ifa")
      param <- "tau"
      scheme <- "horst"
      penalty <- c(1, 1)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair", J = 2)
    },
    "ra" = {
      check_nblocks(blocks, "ra")
      param <- "tau"
      ncomp[2] <- ncomp[1]
      scheme <- "horst"
      penalty <- c(1, 0)
      response <- 2
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "response", J = 2)
    },
    "gcca" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "maxvar" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "maxvar-b" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(0, J + 1)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "maxvar-a" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "mfa" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- rep(1, J + 1)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
      scale_block <- "lambda1"
    },
    "mcia" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      comp_orth <- FALSE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
      scale_block <- "inertia"
    },
    "mcoa" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      comp_orth <- FALSE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
      scale_block <- "inertia"
    },
    "cpca-1" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "horst"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "cpca-2" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- "factorial"
      penalty <- c(rep(1, J), 0)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "cpca-4" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "hpca" = {
      param <- "tau"
      ncomp <- rep(max(ncomp), J + 1)
      scheme <- function(x) x^4
      penalty <- c(rep(1, J), 0)
      response <- NULL
      comp_orth <- TRUE
      superblock <- TRUE
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
    },
    "maxbet-b" = {
      param <- "tau"
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- FALSE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "all")
    },
    "maxbet" = {
      param <- "tau"
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- FALSE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "all")
    },
    "maxdiff-b" = {
      param <- "tau"
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- FALSE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "maxdiff" = {
      param <- "tau"
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- FALSE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "sabscor" = {
      param <- "tau"
      scheme <- "centroid"
      penalty <- rep(0, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "ssqcor" = {
      param <- "tau"
      scheme <- "factorial"
      penalty <- rep(0, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "ssqcov-1" = {
      param <- "tau"
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "all")
    },
    "ssqcov-2" = {
      param <- "tau"
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "ssqcov" = {
      param <- "tau"
      scheme <- "factorial"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "sumcor" = {
      param <- "tau"
      scheme <- "horst"
      penalty <- rep(0, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "sumcov-1" = {
      param <- "tau"
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "all")
    },
    "sumcov-2" = {
      param <- "tau"
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "sumcov" = {
      param <- "tau"
      scheme <- "horst"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
    },
    "sabscov-1" = {
      param <- "tau"
      scheme <- "centroid"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "all")
    },
    "sabscov-2" = {
      param <- "tau"
      scheme <- "centroid"
      penalty <- rep(1, J)
      response <- NULL
      comp_orth <- TRUE
      superblock <- FALSE
      connection <- connection_matrix(blocks, type = "pair")
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
    if (length(modified_parameters) > 0) {
      message(
        "Choice of method '", method, "' overwrote parameters '",
        paste(modified_parameters, collapse = "', '"), "'."
      )
    }
  }

  if (method %in% generic_methods()) {
    scheme <- check_scheme(scheme)
    if (any(sparsity != 1) & (method == "rgcca")) {
      param <- "sparsity"
      method <- "sgcca"
      penalty <- sparsity
    }
    if (!is.null(response)) {
      check_blockx("response", response, blocks)
      ncomp[response] <- max(ncomp[-response])
      superblock <- FALSE
      if (!is.null(connection)) {
        connection <- check_connection(connection, blocks)
      } else {
        connection <- connection_matrix(
          blocks, type = "response", response = response
        )
      }
    }
    penalty <- check_penalty(
      penalty, blocks, method, superblock = superblock, ncomp = max(ncomp)
    )
    if (superblock) {
      ncomp <- rep(max(ncomp), J + 1)
      connection <- connection_matrix(blocks, type = "response", J = J + 1)
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
      if (is.null(connection)) {
        connection <- connection_matrix(blocks, type = "pair")
      } else {
        connection <- check_connection(connection, blocks)
      }
    }
  }
  ncomp <- check_ncomp(
    ncomp, blocks,
    superblock = superblock, response = response
  )

  rgcca_args[[param]] <- penalty

  ### FIX HERE #### -> netsgcca needs other parameters
  if (method == "netsgcca") {
    param_list <- list(lambda = lambda, graph_laplacians = graph_laplacians)
  } else {
    param_list <- list()
  }
  ### end ###

  rgcca_args <- modifyList(rgcca_args, list(
    ncomp = ncomp,
    scheme = scheme,
    method = method,
    response = response,
    comp_orth = comp_orth,
    connection = connection,
    superblock = superblock,
    scale_block = scale_block
  ), keep.null = TRUE)
  return(list(
    rgcca_args = rgcca_args,
    opt = list(
      # gcca = gcca,
      # supplementary_parameters = param_list,
      param = param
    )
  ))
}
