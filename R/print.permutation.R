#' Print a fitted rgcca_permutation object
#'
#' Print a fitted rgcca_permutation object
#' @param x A fitted rgcca_permutation object (see
#' \code{\link[RGCCA]{rgcca_permutation}})
#' @param ... additional print parameters
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

  penalties <- round(x$penalties, 3)
  rownames(penalties) <- seq_len(NROW(penalties))
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", x$par_type, ") used: "), fill = TRUE)
  print(penalties, quote = FALSE, ...)
  cat("\n")

  if (length(x$call$blocks) > 5) {
    combinations <- paste("Tuning parameter set ",
      sep = "",
      seq_along(x$pvals)
    )
  } else {
    combinations <- apply(
      format(x$penalties, digits = 2), 1, paste0,
      collapse = "/"
    )
  }

  tab <- cbind(
    combinations,
    format(cbind(x$crit, x$means, x$sds, x$zstat, x$pvals), digits = 3)
  )
  dimnames(tab) <- list(
    seq_len(NROW(x$penalties)),
    c(
      "Tuning parameters", "Criterion", "Permuted criterion",
      "sd", "zstat", "p-value"
    )
  )
  print(tab, quote = FALSE, ...)

  best <- which(apply(
    x$penalties, 1, function(z) identical(z, x$bestpenalties)
  ))
  cat(paste0(
    "\nThe best combination is: ",
    paste(format(x$bestpenalties, digits = 3), collapse = ", "),
    " for a z score of ", format(x$zstat[best], digits = 3),
    " and a p-value of ", format(x$pvals[best], digits = 3),
    ".\n"
  ))
}
