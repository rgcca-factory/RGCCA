#' Print a fitted object from the RGCCA package
#'
#' @description
#' `print.rgcca()` prints a fitted RGCCA object. The method and number of
#' components are displayed.
#'
#' `print.rgcca_cv()` prints a rgcca_cv object. The type of validation,
#' the number of tried parameter sets, the type of task, and the model used
#' are displayed.
#'
#' `print.rgcca_permutation()` prints a rgcca_permutation object.
#' The number of permutations and tried parameter sets are displayed.
#'
#' `print.rgcca_bootstrap()` prints a rgcca_bootstrap object.
#' The number of boostrap samples used for fitting is displayed.
#'
#' `print.rgcca_stability()` prints a rgcca_stability object.
#' The number of boostrap samples used for fitting is displayed.
#'
#' @param x An object to be printed
#' (output of functions \code{\link{rgcca}},
#' \code{\link{rgcca_cv}}, \code{\link{rgcca_permutation}},
#' \code{\link{rgcca_bootstrap}}, or \code{\link{rgcca_stability}}).
#' @param ... Further arguments passed to other methods.
#' @return none
#'
#' @export
#' @examples
#' ## Printing of an rgcca object
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
#' print(res)
#'
#' ## Printing of an rgcca_cv object
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca", par_type = "tau",
#'   par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1,
#'   verbose = TRUE
#' )
#' print(res)
#'
#' ## Printing of an rgcca_permutation object
#' perm.out <- rgcca_permutation(blocks,
#'   par_type = "tau",
#'   n_perms = 5, n_cores = 1,
#'   verbose = TRUE
#' )
#' print(perm.out)
#'
#' ## Printing of an rgcca_bootstrap object
#' fit.rgcca <- rgcca(blocks, ncomp = c(2, 1, 2))
#' boot.out <- rgcca_bootstrap(fit.rgcca, n_boot = 20, n_cores = 2,
#'                             verbose = TRUE)
#' print(boot.out)
#'
#' ## Printing of an rgcca_stability object
#' fit.sgcca <- rgcca(blocks, sparsity = c(.8, .9, .6))
#' res <- rgcca_stability(fit.sgcca, n_boot = 10, verbose = TRUE)
#' print(res)
#' @export
#' @rdname print
#' @order 1
print.rgcca <- function(x, ...) {
  stopifnot(is(x, "rgcca"))
  cat(
    "Fitted", toupper(x$call$method), "model. \n"
  )
  if (max(x$call$ncomp) == 1) {
    cat(
      "The algorithm converged to a stationnary point after",
      length(x$crit) - 1, "iterations."
    )
  } else {
    cat("The algorithm converged to a stationnary point:")
    for (k in seq_len(max(x$call$ncomp))) {
      cat("\n\t")
      cat(
        "- After ", length(x$crit[[k]]) - 1,
        " iterations for component ", k, ".",
        sep = ""
      )
    }
  }
  cat("\n")
}
