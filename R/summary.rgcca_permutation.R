#' @export
#' @rdname summary
#' @order 3
summary.rgcca_permutation <- function(object, ...) {
  stopifnot(is(object, "rgcca_permutation"))

  ### Print parameters of the function
  print_call(object$call)
 
  params <- object$params
  rownames(params) <- seq_len(NROW(params))
params_clean <- as.data.frame(params, stringsAsFactors = FALSE)

# 2. Loop through each column and collapse any matrix or list data into 'value1/value2'
for (col_name in names(params_clean)) {
  col_data <- params_clean[[col_name]]
  
  if (is.matrix(col_data)) {
    params_clean[[col_name]] <- apply(col_data, 1, paste, collapse = "/")
  } else if (is.list(col_data)) {
    params_clean[[col_name]] <- sapply(col_data, function(x) paste(x, collapse = "/"))
  } else if (is.data.frame(col_data)) {
    params_clean[[col_name]] <- apply(as.matrix(col_data), 1, paste, collapse = "/")
  }
}
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", object$par_type, ") used: "), fill = TRUE)
  print(params_clean, quote = FALSE, ...)
  cat("\n")

  tab <- format(object$stats, digits = 3)
  colnames(tab) <- c(
    "Tuning parameters", "Criterion",
    "Permuted criterion", "sd", "zstat", "p-value"
  )
  print(tab, quote = FALSE, ...)

  best <- which(apply(
    object$params, 1, function(z) identical(z, object$best_params)
  ))
  cat(strwrap(paste0(
    "\nThe best combination is: ",
    object$stats$combinations[best],
    " for a z score of ", format(object$stats$zstat[best], digits = 3),
    " and a p-value of ", format(object$stats$pval[best], digits = 3)
  ), getOption("width")))
  cat("\n")
}
