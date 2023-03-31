#' Print a rgcca_permutation object
#'
#' Print a fitted rgcca_permutation object. Parameters of the
#' analysis, tuning parameters and statistics for each set of
#' parameters are displayed.
#' @param x A fitted rgcca_permutation object (see
#' \code{\link[RGCCA]{rgcca_permutation}})
#' @param ... Other parameters used in print (for the displaying of matrices)
#' @return none
#' @export
#' @examples
#' data(Russett)
#' A <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#'
#' perm.out <- rgcca_permutation(A,
#'   par_type = "tau",
#'   n_perms = 5, n_cores = 1
#' )
#' print(perm.out)
print.permutation <- function(x, ...) {
  ### Print parameters of the function
  print_call(x$call)

  params <- round(x$params, 3)
  rownames(params) <- seq_len(NROW(params))
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", x$par_type, ") used: "), fill = TRUE)
  print(params, quote = FALSE, ...)
  cat("\n")

  tab <- format(x$stats, digits = 3)
  colnames(tab) <- c(
    "Tuning parameters", "Criterion",
    "Permuted criterion", "sd", "zstat", "p-value"
  )
  print(tab, quote = FALSE, ...)

  best <- which(apply(
    x$params, 1, function(z) identical(z, x$best_params)
  ))
  cat(strwrap(paste0(
    "\nThe best combination is: ",
    x$stats$combinations[best],
    " for a z score of ", format(x$stats$zstat[best], digits = 3),
    " and a p-value of ", format(x$stats$pval[best], digits = 3),
    ".\n"
  ), getOption("width")))
}
