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
#' @inheritParams rgcca_bootstrap
#' @inheritParams rgcca
#' @param par_type A character string indicating the parameters to tune between
#' "sparsity" and "tau".
#' @param par_length A numeric value indicating the number of sets of parameters
#' to be tested (if par_value = NULL).
#' @param par_value Sets of parameters to consider during the permutation
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
#' @param n_perms The number of permutations for each set of constraints
#' (default is 20).
#' @param verbose A logical value indicating if the progress of the
#' permutation procedure is reported.
#' @return  \item{opt}{A list indicating some options of the RGCCA model used
#' during the permutation.}
#' @return \item{call}{A list containing the input parameters}
#' @return \item{par_type}{A character string indicating the parameters to tune
#' between "sparsity" and "tau".}
#' @return \item{n_perms}{The number of permutations for each set of candidate
#' tuning paramaters (default
#' is 20)}
#' @return \item{best_params}{The set of tuning parameters that yields the
#' highest Z-statistic}
#' @return \item{permcrit}{A matrix of permuted S/RGCCA criteria. The ith row of
#' permcrit contains the n_perms values of S/RGCCA permuted criteria
#' obtained for each set of tuning parameters.}
#' @return \item{params}{A matrix giving, the set of candidate tuning
#' paramaters considered during the permutation process (tau or sparsity).}
#' @return \item{stats}{A data.frame containing in column: the set of candidate
#' parameters, the corresponding non permuted criterion, mean and standard
#' deviation of permuted criteria, Z-statistic and p-value.}
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
#' fit$best_params
#'
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
#' perm.out$best_params
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
#' perm.out$best_params
#'
#' # when par_value is a scalar, the same maximum value is applied
#' # for each block
#'
#' perm.out <- rgcca_permutation(blocks,
#'   connection = C,
#'   par_type = "sparsity",
#'   par_value = 0.8, par_length = 5,
#'   n_perms = 10, n_cores = 1
#' )
#'
#' perm.out$params
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
#' data("ge_cgh_locIGR", package = "gliomaData")
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
                              connection = NULL,
                              scheme = "factorial",
                              ncomp = 1,
                              tau = 1,
                              sparsity = 1,
                              init = "svd", bias = TRUE, tol = 1e-8,
                              response = NULL, superblock = FALSE,
                              NA_method = "nipals", rgcca_res = NULL,
                              verbose = TRUE, n_iter_max = 1000,
                              comp_orth = TRUE) {
  ### Try to retrieve parameters from a rgcca object
  rgcca_args <- as.list(environment())
  tmp <- get_rgcca_args(blocks, rgcca_args)
  opt <- tmp$opt
  rgcca_args <- tmp$rgcca_args

  ### Check parameters
  check_integer("n_perms", n_perms)
  check_integer("par_length", n_perms)
  match.arg(par_type, c("tau", "sparsity", "ncomp"))
  if (length(rgcca_args$blocks) == 1) {
    stop_rgcca(
      "wrong number of blocks. Permutation requires more than ",
      "one block."
    )
  }

  ### Prepare parameters for line search
  if (
    rgcca_args$method %in% sparse_methods() && (par_type == "tau")
  ) {
    par_type <- "sparsity"
  } else if (par_type == "sparsity") {
    rgcca_args$method <- "sgcca"
    opt$gcca <- sgcca
    opt$param <- "sparsity"
  }

  param <- set_parameter_grid(
    par_type, par_length, par_value, rgcca_args$blocks,
    rgcca_args[[par_type]], rgcca_args$response,
    rgcca_args$superblock,  opt$disjunction
  )

  # Generate a warning if tau has not been fully specified for a block that
  # has more columns than samples and remove tau = 0 configuration
  n <- NROW(rgcca_args$blocks[[1]])
  overfitting_risk <- (param$par_type == "tau") && is.null(dim(par_value)) &&
    any(vapply(
      seq_along(rgcca_args$blocks),
      function(j) NCOL(rgcca_args$blocks[[j]]) > n,
      FUN.VALUE = logical(1L)
    ))
  if (overfitting_risk) {
    param$par_value <- param$par_value[-nrow(param$par_value), ]
    warning(
      "overfitting risk. A block has more columns than rows, so the ",
      "configuration with tau = 0 has been removed."
    )
  }

  ### Create folds
  v_inds <- lapply(seq_len(n_perms), function(i) {
    lapply(rgcca_args$blocks, function(x) {
      sample(seq_len(NROW(x)))
    })
  })

  ### Start line search
  # For every set of parameter, RGCCA is run once on the non permuted blocks
  # and then n_perms on permuted blocks.
  idx <- seq(NROW(param$par_value) * (n_perms + 1))
  W <- par_pblapply(idx, function(n) {
    i <- (n - 1) %/% (n_perms + 1) + 1
    j <- (n - 1) %% (n_perms + 1)
    perm <- (n - 1) %% (n_perms + 1) != 0
    rgcca_permutation_k(
      rgcca_args,
      inds = v_inds[[j]],
      par_type = param$par_type,
      par_value = param$par_value[i, ],
      perm = perm
    )
  }, n_cores = n_cores, verbose = verbose)

  W <- do.call(rbind, W)

  ### Format output
  par_colnames <- names(rgcca_args$blocks)
  if (ncol(param$par_value) > length(rgcca_args$blocks)) {
    par_colnames <- c(par_colnames, "superblock")
  }
  rownames(param$par_value) <- seq_len(NROW(param$par_value))
  colnames(param$par_value) <- par_colnames

  idx_perm <- (idx - 1) %% (n_perms + 1) != 0
  crit <- W[!idx_perm]
  permcrit <- matrix(W[idx_perm],
    nrow = nrow(param$par_value),
    ncol = n_perms, byrow = TRUE
  )

  # Compute statistics
  pvals <- vapply(
    seq_len(NROW(param$par_value)),
    function(k) mean(permcrit[k, ] >= crit[k]),
    FUN.VALUE = double(1)
  )
  zstat <- vapply(
    seq_len(NROW(param$par_value)), function(k) {
      z <- (crit[k] - mean(permcrit[k, ])) / (sd(permcrit[k, ]))
      if (is.na(z) || z == "Inf") {
        z <- 0
      }
      return(z)
    },
    FUN.VALUE = double(1)
  )
  combinations <- format_combinations(param$par_value)

  stats <- data.frame(
    combinations = combinations,
    crit = crit,
    mean = apply(permcrit, 1, mean, na.rm = TRUE),
    sd = apply(permcrit, 1, sd, na.rm = TRUE),
    zstat = zstat,
    pval = pvals
  )

  structure(list(
    opt = opt, call = rgcca_args, par_type = par_type,
    n_perms = n_perms, best_params = param$par_value[which.max(zstat), ],
    permcrit = permcrit, params = param$par_value, stats = stats
  ), class = "permutation")
}
