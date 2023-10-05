#' Summary of a fitted object from the RGCCA package
#'
#' @description
#' `summary.rgcca()` summarizes a fitted RGCCA object.
#' Some information about the
#' model are displayed like model parameters or criterion.
#'
#' `summary.rgcca_cv()` summarizes a fitted rgcca_cv object. Parameters of the
#' analysis, tuning parameters and statistics for each set of
#' parameters are displayed.
#'
#' `summary.rgcca_permutation()` summarizes a fitted rgcca_permutation object.
#' Parameters of the analysis, tuning parameters and statistics for each set of
#' parameters are displayed.
#'
#' `summary.rgcca_bootstrap()` summarizes a fitted rgcca_bootstrap object.
#' Parameters of the analysis and bootstrap statistics are displayed.
#'
#' `summary.rgcca_stability()` calls `summary.rgcca()` on the fitted RGCCA model
#' returned by `rgcca_stability()`.
#'
#' @inheritParams plot.rgcca
#' @param object An object to be summarized
#' (output of functions \code{\link{rgcca}},
#' \code{\link{rgcca_cv}}, \code{\link{rgcca_permutation}},
#' \code{\link{rgcca_bootstrap}}, or \code{\link{rgcca_stability}}).
#' @param type A character string indicating the type of the summarized object
#' (see details).
#' @param block A numeric corresponding to the block(s) to summarize.
#' @param ... Further arguments passed to other methods
#' (for the displaying of matrices).
#' @return none
#' @details
#' Argument type can take two values in `summary.cval`: \itemize{
#' \item "sd" (default): mean values of the cross-validated scores are reported,
#' as well as means plus or minus standard deviations.
#' \item "quantiles": median values, 25\% and 75\% quantiles of the
#' cross-validated scores are reported.
#' }
#'
#' Argument type can take two values in `summary.bootstrap`: \itemize{
#' \item "weights" (default): statistics about the block-weight vectors
#' are reported.
#' \item "loadings": statistics about the block-loading vectors are reported.
#' }
#'
#' @export
#' @examples
#' ## Summary of an rgcca object
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:8]
#' )
#' C <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' res <- rgcca(blocks,
#'   connection = C, ncomp = rep(2, 3), tau = c(1, 1, 1),
#'   scheme = "factorial", scale = TRUE, verbose = FALSE
#' )
#' summary(res)
#'
#' ## Summary of an rgcca_cv object
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca", par_type = "tau",
#'   par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1,
#'   verbose = TRUE
#' )
#' summary(res)
#'
#' ## Summary of an rgcca_permutation object
#' perm.out <- rgcca_permutation(blocks,
#'   par_type = "tau",
#'   n_perms = 5, n_cores = 1,
#'   verbose = TRUE
#' )
#' summary(perm.out)
#'
#' ## Summary of an rgcca_bootstrap object
#' fit.rgcca <- rgcca(blocks, ncomp = c(2, 1, 2))
#' boot.out <- rgcca_bootstrap(fit.rgcca, n_boot = 20, n_cores = 2,
#'                             verbose = TRUE)
#' summary(boot.out)
#'
#' ## Summary of an rgcca_stability object
#' fit.sgcca <- rgcca(blocks, sparsity = c(.8, .9, .6))
#' res <- rgcca_stability(fit.sgcca, n_boot = 10, verbose = TRUE)
#' summary(res)
#' @rdname summary
#' @order 1
summary.rgcca <- function(object, ...) {
  ### Print parameters of the function
  print_call(object$call)

  ### Print criterion
  if (is.list(object$crit)) {
    crit <- Reduce("+", lapply(object$crit, function(t) {
      return(t[length(t)])
    }))
  } else {
    crit <- object$crit[length(object$crit)]
  }
  cat("Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = ",
    sep = "",
    paste(round(crit, 4), sep = "", " "), fill = TRUE
  )

  ### Print regularization parameter or the number of selected variables
  cat("\n")
  if (!tolower(object$call$method) %in% sparse_methods()) {
    param <- "regularization"
    if (!is.matrix(object$call$tau)) {
      for (i in seq_len(NCOL(object$call$connection))) {
        tau <- object$call$tau[i]
        cat("The", param, "parameter used for", names(object$blocks)[i],
          "is:", round(tau, 4),
          fill = TRUE
        )
      }
    } else {
      cat("The", param, "parameters used are: \n")
      print(round(object$call$tau, 4), ...)
    }
  } else {
    response <- ifelse(
      object$opt$disjunction, object$call$response, length(object$blocks) + 1
    )
    nb_selected_var <- lapply(
      object$a[-response],
      function(a) apply(a, 2, function(l) sum(l != 0))
    )
    param <- "sparsity"
    if (!is.matrix(object$call$sparsity)) {
      for (i in seq_len(NCOL(object$call$connection))[-response]) {
        sparsity <- object$call$sparsity[i]

        cat("The", param, "parameter used for", names(object$blocks)[i], "is:",
          sparsity, "(with", paste(nb_selected_var[[i]], collapse = ", "),
          "variables selected)",
          fill = TRUE
        )
      }
    } else {
      cat("The", param, "parameters used are: \n")
      print(round(object$call$sparsity[, -response], 4), ...)
      cat("The number of selected variables are: \n")
      print(do.call(cbind, nb_selected_var))
    }
    if (object$opt$disjunction) {
      cat("The regularization parameter used for",
          names(object$blocks)[response], "is:", 0,
          fill = TRUE
      )
    }
  }
  cat("\n")
}
