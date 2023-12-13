#' @export
#' @rdname summary
#' @order 3
summary.rgcca_permutation <- function(object, ...) {
  stopifnot(is(object, "rgcca_permutation"))

  ### Print parameters of the function
  print_call(object$call)

  params <- round(object$params, 3)
  rownames(params) <- seq_len(NROW(params))
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", object$par_type, ") used: "), fill = TRUE)
  print(params, quote = FALSE, ...)
  cat("\n")

  tab <- format(object$stats, digits = 3)
  colnames(tab) <- c(
    "Tuning parameters", "Criterion",
    "Permuted criterion", "sd", "zstat", "p-value"
  )
  print(tab, quote = FALSE, ...)

  best <- which(apply(
    object$params, 1, function(z) identical(z, drop(object$best_params))
  ))
  cat(strwrap(paste0(
    "\nThe best combination is: ",
    object$stats$combinations[best],
    " for a z score of ", format(object$stats$zstat[best], digits = 3),
    " and a p-value of ", format(object$stats$pval[best], digits = 3)
  ), getOption("width")))
  cat("\n")
}
