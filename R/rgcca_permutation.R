#' Tune the S/RGCCA hyper-parameters by permutation
#'
#' This function can be used to automatically select the hyper-parameters
#' (amount of sparsity for sgcca or shrinkage parameters for RGCCA).
#' A permutation based strategy very similar to the one proposed in
#' (Witten et al, 2009) is implemented.
#'
#' @details
#' The tuning parameters are selected using the permutation scheme proposed in
#' (Witten et al, 2009). For each candidate tuning parameter value, the
#' following is performed:
#'
#' (1) Repeat the following n_perms times (for n_perms large): \cr
#'    \verb{    }(a) The samples in \eqn{X_1},..., \eqn{X_J} are randomly
#'    permuted blocks: \eqn{X_1^*},..., \eqn{X_J^*}. \cr
#'    \verb{    }(b) S/RGCCA is run on the permuted data sets \eqn{X_1^*},...,
#'       \eqn{X_J^*} to get canonical variates \eqn{a_1^*},..., \eqn{a_J^*}.\cr
#'    \verb{    }(c) Record t* = sum_(j,k) c_jk g(Cov(X_j^*a_j^*, X_k^*a_k^*).
#'
#' (2) Sparse CCA is run on the blocks \eqn{X_1},..., \eqn{X_J} to obtain
#'     canonical variates \eqn{a_1},..., \eqn{a_J}.
#'
#' (3) Record t = sum_(j,k) c_jk g(Cov(X_ja_j, X_ka_k).
#'
#' (4) The resulting p-value is given by $mean(t* > t)$; that is, the fraction
#' of t* that exceed the value of t obtained from the real data.
#'
#' Then, choose the tuning parameter values that gives the smallest value in
#' Step 4.
#'
#' This function only selects tuning parameters for the first deflation stage
#' of S/RGCCA. By default, this function performs a one-dimensional
#' search in tuning parameter space.
#'
#' @inheritParams bootstrap
#' @inheritParams rgcca
#' @inheritParams plot2D
#' @param par_type A character string indicating the parameters to tune between
#' "sparsity" and "tau".
#' @param par_length A numeric value indicating the number of sets of penalties
#' to be tested (if par_value = NULL).
#' @param par_value Sets of penalties to consider during the permutation
#' process.
#' If par_value = NULL, it takes 10 sets between min values (0 for RGCCA and
#' 1/sqrt(ncol(Xj)) for SGCCA) and 1. Otherwise, it could be either (i) A matrix
#' of dimension IxJ (where I the number of combinations to be tested and J the
#' number of blocks), or (ii) a vector of length J length specifying the maximal
#' values to consider for each block. In that case, par_length combinations
#' are tested from min values to the maximal values specified by this vector.
#' (iii) a numerical value giving the same maximal value to be considered for
#' each block. In that case par_length combinations are tested from min values
#' to this single maximal value.
#' @param n_perms Number of permutations for each set of constraints (default
#' is 20).
#' @param verbose Logical value indicating if the progress of the
#' permutation procedure is reported.
#' @return \item{zstat}{A vector of Z-statistics, one zstat per set of tuning
#' parameters.}
#' @return \item{bestpenalties}{The set of tuning parameters that yields the
#' highest Z-statistics}
#' @return \item{permcrit}{Matrix of permuted S/RGCCA criteria. The ith row of
#' permcrit contains the n_perms values of S/RGCCA permuted criteria
#' obtained for each set of tuning parameters.}
#' @return \item{means}{A vector that contains, for each set of tuning
#' parameters, the mean of the permuted R/SGCCA criteria}
#' @return \item{sds}{A vector that contains, for each set of tuning
#' parameters, the standard deviation of the permuted R/SGCCA criteria}
#' @return \item{crit}{A vector that contains, for each set of tuning
#' parameters, the value of the R/SGCCA criteria obtained from the original
#' data.}
#' @return \item{pval}{Vector of p-values, one per set of tuning parameters.}
#' @return \item{penalties}{Matrix giving, the set of tuning paramaters
#' considered during the permutation process (tau or sparsity).}
#' @references Witten, D. M., Tibshirani, R., & Hastie, T. (2009). A penalized
#' matrix decomposition, with applications to sparse principal components and
#' canonical correlation analysis. Biostatistics, 10(3), 515-534.
#' @examples
#' ####################################
#' # Permutation based strategy for   #
#' # determining the best shrinkage   #
#' # parameters (par_type = "tau")    #
#' ####################################
#'
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#'
#' C <- matrix(c(
#'   0, 0, 1,
#'   0, 0, 1,
#'   1, 1, 0
#' ), 3, 3)
#'
#' # default value: 10 vectors from rep(0, length(blocks))
#' # to rep(1, length(blocks)), uniformly distributed.
#'
#' fit <- rgcca_permutation(blocks,
#'   connection = C,
#'   par_type = "tau",
#'   par_length = 10, n_perms = 2,
#'   n_cores = 1
#' )
#'
#' print(fit)
#' plot(fit)
#' fit$bestpenalties
#' \dontrun{
#' # It is possible to define explicitly K combinations of shrinkage
#' # parameters to be tested and in that case a matrix of dimension KxJ is
#' # required. Each row of this matrix corresponds to one specific set of
#' # shrinkage parameters.
#' par_value <- matrix(c(
#'   0, 0, 0,
#'   1, 1, 0,
#'   0.5, 0.5, 0.5,
#'   sapply(blocks, RGCCA:::tau.estimate),
#'   1, 1, 1
#' ), 5, 3, byrow = TRUE)
#'
#'
#' perm.out <- rgcca_permutation(blocks,
#'   connection = C,
#'   par_type = "tau",
#'   par_value = par_value,
#'   n_perms = 5, n_cores = 1
#' )
#'
#' print(perm.out)
#' plot(perm.out)
#'
#' # with superblock
#'
#' perm.out <- rgcca_permutation(blocks,
#'   par_type = "tau",
#'   superblock = TRUE,
#'   scale = TRUE, scale_block = FALSE,
#'   n_perms = 5, n_cores = 1
#' )
#'
#' print(perm.out)
#' plot(perm.out)
#'
#' # used a fitted rgcca_permutation object as input of the rgcca function
#' fit.rgcca <- rgcca(perm.out)
#' fit.rgcca$call$tau
#' fit.rgcca$call$scale
#' fit.rgcca$call$scale_block
#'
#' ######################################
#' # Permutation based strategy for     #
#' # determining the best sparsity      #
#' # parameters (par_type = "sparsity") #
#' ######################################
#'
#' # defaut value: 10 vectors from minimum values
#' # (1/sqrt(ncol(X1)), ..., 1/sqrt(ncol(XJ))
#' # to rep(1, J), uniformly distributed.
#'
#' perm.out <- rgcca_permutation(blocks,
#'   par_type = "sparsity",
#'   n_perms = 50, n_cores = 1
#' )
#'
#' print(perm.out)
#' plot(perm.out)
#' perm.out$bestpenalties
#'
#' # when par_value is a vector of length J. Each element of the vector
#' # indicates the maximum value of sparsity to be considered for each block.
#' # par_length (default value = 10) vectors from minimum values
#' # (1/sqrt(ncol(X1)), ..., 1/sqrt(ncol(XJ)) to maximum values, uniformly
#' # distributed, are then considered.
#'
#' perm.out <- rgcca_permutation(blocks,
#'   connection = C,
#'   par_type = "sparsity",
#'   par_value = c(0.6, 0.75, 0.5),
#'   par_length = 7, n_perms = 20,
#'   n_cores = 1, tol = 1e-3
#' )
#'
#' print(perm.out)
#' plot(perm.out)
#' perm.out$bestpenalties
#'
#' # when par_value is a scalar, the same maximum value is applied
#' # for each block
#'
#' perm.outt <- rgcca_permutation(blocks,
#'   connection = C,
#'   par_type = "sparsity",
#'   par_value = 0.8, par_length = 5,
#'   n_perms = 10, n_cores = 1
#' )
#'
#' perm.out$penalties
#'
#' ######################################
#' # speed up the permutation procedure #
#' ######################################
#'
#' # The rgcca_permutation function can be quite time-consuming. Since
#' # approximate estimates of the block weight vectors are acceptable in this
#' # case, it is possible to reduce the value of the tolerance (tol argument)
#' # of the RGCCA algorithm to speed up the permutation procedure.
#' #
#' require(gliomaData)
#' data(ge_cgh_locIGR)
#' A <- ge_cgh_locIGR$multiblocks
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' A[[3]] <- A[[3]][, -3]
#' C <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#'
#' # check dimensions of the blocks
#' sapply(A, dim)
#'
#' par_value <- matrix(c(
#'   seq(0.1, 1, by = 0.1),
#'   seq(0.1, 1, by = 0.1),
#'   rep(0, 10)
#' ), 10, 3, byrow = FALSE)
#'
#' fit <- rgcca_permutation(A,
#'   connection = C,
#'   par_type = "tau",
#'   par_value = par_value,
#'   par_length = 10,
#'   n_perms = 10, n_cores = 1, tol = 1e-2
#' )
#' print(fit)
#' plot(fit)
#' }
#'
#' @export
rgcca_permutation <- function(blocks, par_type = "tau", par_value = NULL,
                              par_length = 10, n_perms = 20,
                              n_cores = 1,
                              quiet = TRUE, scale = TRUE, scale_block = TRUE,
                              method = "rgcca",
                              connection = 1 - diag(length(blocks)),
                              scheme = "factorial",
                              ncomp = rep(1, length(blocks)),
                              tau = rep(1, length(blocks)),
                              sparsity = rep(1, length(blocks)),
                              init = "svd", bias = TRUE, tol = 1e-8,
                              response = NULL, superblock = FALSE,
                              NA_method = "nipals", rgcca_res = NULL,
                              verbose = TRUE) {
  ### Try to retrieve parameters from a rgcca object
  if (!missing(blocks) & class(blocks) == "rgcca") {
    rgcca_res <- blocks
  }
  if (class(rgcca_res) == "rgcca") {
    stopifnot(is(rgcca_res, "rgcca"))
    message("All parameters were imported from a fitted rgcca object.")
    scale_block <- rgcca_res$call$scale_block
    scale <- rgcca_res$call$scale
    scheme <- rgcca_res$call$scheme
    response <- rgcca_res$call$response
    tol <- rgcca_res$call$tol
    NA_method <- rgcca_res$call$NA_method
    init <- rgcca_res$call$init
    bias <- rgcca_res$call$bias
    blocks <- rgcca_res$call$raw
    superblock <- rgcca_res$call$superblock
    connection <- rgcca_res$call$connection
    tau <- rgcca_res$call$tau
    ncomp <- rgcca_res$call$ncomp
    sparsity <- rgcca_res$call$sparsity
    method <- rgcca_res$call$method
    superblock <- rgcca_res$call$superblock
  }

  ### Check parameters
  check_integer("n_perms", n_perms)
  check_integer("par_length", n_perms)
  match.arg(par_type, c("tau", "sparsity", "ncomp"))
  if (length(blocks) == 1) {
    stop_rgcca(
      "wrong number of blocks. Permutation requires more than ",
      "one block."
    )
  }

  ### Prepare parameters for line search
  if (method %in% c("sgcca", "spca", "spls") && (par_type == "tau")) {
    par_type <- "sparsity"
  } else if (par_type == "sparsity") method <- "sgcca"

  call <- list(
    method = method, par_type = par_type, par_value = par_value,
    n_perms = n_perms, quiet = quiet, connection = connection,
    NA_method = NA_method, tol = tol, scheme = scheme,
    scale = scale, scale_block = scale_block,
    superblock = superblock, blocks = blocks, ncomp = ncomp,
    tau = tau, sparsity = sparsity, response = response
  )

  param <- set_parameter_grid(
    par_type, par_length, par_value, blocks, response,
    superblock
  )

  ### Start line search
  # For every set of parameter, RGCCA is run once on the non permuted blocks
  # and then n_perms on permuted blocks.
  idx <- seq(NROW(param$par_value) * (n_perms + 1))
  W <- unlist(par_pblapply(idx, function(n) {
    i <- (n - 1) %/% (n_perms + 1) + 1
    perm <- (n - 1) %% (n_perms + 1) != 0
    rgcca_permutation_k(
      blocks = blocks,
      par_type = param$par_type,
      par_value = param$par_value[i, ],
      perm = perm,
      method = method,
      quiet = quiet,
      superblock = superblock,
      scheme = scheme,
      tol = tol,
      scale = scale,
      scale_block = scale_block,
      connection = connection,
      NA_method = NA_method,
      bias = bias,
      init = init,
      ncomp = ncomp,
      tau = tau,
      sparsity = sparsity
    )
  }, n_cores = n_cores, verbose = verbose))

  ### Format output
  par_colnames <- names(blocks)
  if (ncol(param$par_value) > length(blocks)) {
    par_colnames <- c(par_colnames, "superblock")
  }
  rownames(param$par_value) <- seq(NROW(param$par_value))
  colnames(param$par_value) <- par_colnames

  idx_perm <- (idx - 1) %% (n_perms + 1) != 0
  crits <- W[!idx_perm]
  permcrit <- matrix(W[idx_perm],
    nrow = nrow(param$par_value),
    ncol = n_perms, byrow = TRUE
  )
  means <- apply(permcrit, 1, mean, na.rm = TRUE)
  sds <- apply(permcrit, 1, sd, na.rm = TRUE)
  pvals <- vapply(
    seq(NROW(param$par_value)),
    function(k) mean(permcrit[k, ] >= crits[k]),
    FUN.VALUE = double(1)
  )
  zs <- vapply(
    seq(NROW(param$par_value)), function(k) {
      z <- (crits[k] - mean(permcrit[k, ])) / (sd(permcrit[k, ]))
      if (is.na(z) || z == "Inf") {
        z <- 0
      }
      return(z)
    },
    FUN.VALUE = double(1)
  )

  structure(list(
    call = call, zstat = zs,
    bestpenalties = param$par_value[which.max(zs), ],
    permcrit = permcrit, means = means, sds = sds,
    crit = crits, pvals = pvals, penalties = param$par_value
  ), class = "permutation")
}
