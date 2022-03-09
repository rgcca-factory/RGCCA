#' Define the parameters associated with each multi-block component
#' method of the literature.
#'
#' @param method A character string indicating the multi-block component
#' method to consider: rgcca, sgcca, pca, spca, pls, spls, cca,
#' ifa, ra, gcca, maxvar, maxvar-b, maxvar-a, mcoa,cpca-1, cpca-2,
#' cpca-4, hpca, maxbet-b, maxbet, maxdiff-b, maxdiff, maxvar-a,
#' sabscor, ssqcor, ssqcor, ssqcov-1, ssqcov-2, ssqcov, sumcor,
#' sumcov-1, sumcov-2, sumcov, sabscov, sabscov-1, sabscov-2.
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
#'  design explicitely the sheme function (e.g. function(x) x^4) as argument of
#'  rgcca function.  See (Tenenhaus et al, 2017) for details.
#' @param verbose Logical value indicating whether the warnings are displayed.
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
                            verbose = TRUE,
                            quiet = FALSE,
                            response = NULL) {
  substitute <- function(x, env = environment()) {
    x <- strsplit(deparse(base::substitute(x, env)), "\\$")[[1]]
    return(x[length(x)])
  }
  local_env <- new.env()
  local_env$ncomp <- ncomp
  local_env$blocks <- blocks
  local_env$scheme <- scheme
  local_env$penalty <- penalty
  local_env$superblock <- superblock
  local_env$connection <- connection
  local_env$warn_method_value <-
    local_env$warn_method_par <-
    local_env$warn_message_super <- character(0)
  J <- length(local_env$blocks)
  msg_superblock <- "A superblock is considered."
  msg_type <- paste0("By using a ", toupper(method), ", ")

  if (quiet) {
    verbose <- FALSE
  }

  ### SETTINGS ###

  warn_param <- function(x, param) {
    local_env$warn_method_par <- c(
      local_env$warn_method_par, paste(substitute(x, environment()))
    )
    local_env$warn_method_value <- c(
      local_env$warn_method_value, toString(param)
    )
  }

  set_penalty <- function(x) {
    warn_param(local_env$penalty, x)
    return(x)
  }

  set_scheme <- function(x) {
    warn_param(local_env$scheme, x)
    return(x)
  }

  set_connection <- function(x) {
    warn_param(local_env$connection, paste(substitute(x, environment())))
    return(x)
  }

  warn_super <- function(x) {
    if (class(x) %in% c("matrix", "data.frame") &&
      NCOL(x) < (length(local_env$blocks)) &&
      is.null(response)) {
      local_env$warn_msg_super <- c(
        local_env$warn_msg_super, substitute(x, environment())
      )
      return(cbind(x, 1))
    } else if (length(x) < (length(local_env$blocks)) && is.null(response)) {
      local_env$warn_msg_super <- c(
        local_env$warn_msg_super, substitute(x, environment())
      )
      if (substitute(x, environment()) == "ncomp") {
        return(c(x, max(x)))
      } else {
        return(c(x, 1))
      }
    } else {
      return(x)
    }
  }

  set_superblock <- function(verbose = TRUE) {
    local_env$blocks <- c(
      local_env$blocks,
      superblock = list(Reduce(cbind, local_env$blocks))
    )
    local_env$superblock <- TRUE
    local_env$connection <- NULL
    local_env$ncomp <- warn_super(local_env$ncomp)
  }

  set_2block <- function(method) {
    check_nblocks(local_env$blocks, method)

    local_env$scheme <- set_scheme("horst")
    local_env$connection <- set_connection(1 - diag(2))
  }

  ### CHECK TYPES ###

  if (length(grep("[sr]gcca", tolower(method))) == 1) {
    if (local_env$superblock) {
      set_superblock(FALSE)
      local_env$penalty <- warn_super(local_env$penalty)
    } else {
      local_env$superblock <- FALSE
    }
  } else {
    local_env$superblock <- FALSE
  }

  if (length(grep("^s?pca$", tolower(method))) == 1) {
    check_nblocks(local_env$blocks, method)

    local_env$scheme <- set_scheme("horst")
    set_superblock()
    if (tolower(method) == "pca") {
      local_env$penalty <- set_penalty(c(1, 1))
    }
  }

  # 2 Blocks cases
  else if (tolower(method) %in% c("cca", "ra", "ifa", "pls", "spls")) {
    set_2block(method)

    if (tolower(method) == "cca") {
      local_env$penalty <- set_penalty(c(0, 0))
    } else if (tolower(method) %in% c("ifa", "pls")) {
      local_env$penalty <- set_penalty(c(1, 1))
    } else if (tolower(method) == "ra") {
      local_env$penalty <- set_penalty(c(1, 0))
    }
  }

  # Design matrix of 1 values everywhere
  else if (tolower(method) %in% c(
    "sumcor",
    "ssqcor",
    "sabscor",
    "sumcov-1",
    "maxbet",
    "ssqcov-1",
    "maxbet-b",
    "sabscov-1"
  )) {
    local_env$connection <- set_connection(matrix(1, J, J))

    # COR models
    if (tolower(method) %in% c("sumcor", "ssqcor", "sabscor")) {
      local_env$penalty <- set_penalty(rep(0, J))

      switch(tolower(method),
        "sumcor" = {
          local_env$scheme <- set_scheme("horst")
        },
        "ssqcor" = {
          local_env$scheme <- set_scheme("factorial")
        },
        "sabscor" = {
          local_env$scheme <- set_scheme("centroid")
        }
      )
    }

    # COV models
    else if (tolower(method) %in% c(
      "sumcov-1",
      "maxbet",
      "ssqcov-1",
      "maxbet-b",
      "sabscov-1"
    )) {
      local_env$penalty <- set_penalty(rep(1, J))

      if (tolower(method) %in% c("sumcov-1", "maxbet")) {
        local_env$scheme <- set_scheme("horst")
      } else if (tolower(method) %in% c("ssqcov-1", "maxbet-b")) {
        local_env$scheme <- set_scheme("factorial")
      } else if (tolower(method) %in% c("sabscov-1")) {
        local_env$scheme <- set_scheme("centroid")
      }
    }

    # Design matrix with 1 values everywhere except
    # on the diagonal equals to 0
  } else if (tolower(method) %in% c(
    "sumcov",
    "sumcov-2",
    "maxdiff",
    "ssqcov",
    "ssqcov-2",
    "maxdiff-b",
    "sabscov-2"
  )) {
    local_env$connection <- set_connection(1 - diag(J))

    if (tolower(method) %in% c("sumcov", "sumcov-2", "maxdiff")) {
      local_env$scheme <- set_scheme("horst")
      local_env$penalty <- set_penalty(rep(1, J))
    } else if (tolower(method) %in% c("ssqcov", "ssqcov-2", "maxdiff-b")) {
      local_env$scheme <- set_scheme("factorial")
      local_env$penalty <- set_penalty(rep(1, J))
    }
  }

  # Models with a superblock
  else if (tolower(method) %in% c(
    "gcca", "maxvar", "maxvar-b",
    "cpca-1", "cpca-2", "maxvar-a", "mcoa",
    "cpca-4", "hpca"
  )) {
    set_superblock()

    if (tolower(method) %in% c("gcca", "maxvar", "maxvar-b")) {
      local_env$scheme <- set_scheme("factorial")
      local_env$penalty <- set_penalty(rep(0, J + 1))
    } else if (tolower(method) == "cpca-1") {
      local_env$scheme <- function(x) x
      local_env$penalty <- set_penalty(c(rep(1, J), 0))
    } else if (tolower(method) %in% c("maxvar-a", "cpca-2", "mcoa")) {
      local_env$scheme <- set_scheme("factorial")
      local_env$penalty <- set_penalty(c(rep(1, J), 0))
    } else if (tolower(method) %in% c("hpca", "cpca-4")) {
      local_env$scheme <- function(x) x^4
      local_env$penalty <- set_penalty(c(rep(1, J), 0))
    }
  }

  ### WARNINGS ###
  n <- length(local_env$warn_method_par)
  if (verbose & n > 0) {
    set_plural <- function(x = local_env$warn_method_par,
                           y = local_env$warn_method_value,
                           sep = " and ") {
      local_env$warn_method_par <- paste0(x, collapse = sep)
      local_env$warn_method_value <- paste0(y, collapse = sep)
    }

    if (n > 1) {
      grammar <- "s were respectively"
      if (n == 2) {
        set_plural()
      } else {
        warn.method <- c(
          local_env$warn_method_par[n], local_env$warn_method_value[n]
        )
        set_plural(
          local_env$warn_method_par[-n], local_env$warn_method_value[-n], ", "
        )
        set_plural(
          c(local_env$warn_method_par, warn.method[1]),
          c(local_env$warn_method_value, warn.method[2])
        )
      }
    } else {
      grammar <- " was"
    }

    msg <- paste0(
      local_env$warn_method_par, " parameter",
      grammar, " set to ", local_env$warn_method_value
    )

    if (local_env$superblock & tolower(method) != "pca") {
      msg <- paste0(msg, " and ", msg_superblock)
    }

    warning(paste0(msg_type, msg, "."))
  }

  if (verbose & local_env$superblock) {
    if (n < 0) paste0(msg_superblock, msg_superblock)
  }

  if (!quiet & length(local_env$warn_message_super) > 0) {
    if (length(local_env$warn_message_super) > 1) {
      local_env$warn_message_super <- paste(
        local_env$warn_message_super,
        collapse = " and "
      )
      grammar <- "were those"
    } else {
      grammar <- "was the one"
    }
  }

  return(list(
    scheme = local_env$scheme,
    penalty = local_env$penalty,
    ncomp = local_env$ncomp,
    connection = local_env$connection,
    superblock = local_env$superblock
  ))
}
