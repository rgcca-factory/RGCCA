#' @export
#' @rdname print
#' @order 3
print.rgcca_permutation <- function(x, ...) {
  stopifnot(is(x, "rgcca_permutation"))

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
