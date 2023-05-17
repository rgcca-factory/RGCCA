#' @export
#' @rdname print
#' @order 2
print.rgcca_cv <- function(x, type = c("sd", "quantile"), ...) {
  stopifnot(is(x, "rgcca_cv"))
  type <- type[1]
  type <- match.arg(type, c("sd", "quantile"))

  ### Print parameters of the function
  print_call(x$call)

  params <- round(x$params, 3)
  rownames(params) <- seq_len(NROW(params))
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", x$par_type, ") used: "), fill = TRUE)
  print(params, quote = FALSE, ...)
  cat("\n")

  cat(paste0(
    "Validation: ", x$validation,
    ifelse(x$validation == "kfold",
      paste0(
        " with ", x$k, " folds and ",
        x$n_run, " run(s))"
      ), ""
    )
  ), "\n")
  cat(paste("Prediction model:", x$prediction_model, "\n"))
  cat("\n")

  df <- format(x$stats, digits = 3)
  if (type == "quantile") {
    df <- df[, c(1, 4, 5, 6)]
    colnames(df) <- c(
      "Tuning parameters", paste("Median", x$metric), "Q1", "Q3"
    )
  } else {
    df <- df[, c(1, 2, 3)]
    colnames(df) <- c("Tuning parameters", paste("Mean", x$metric), "Sd")
  }
  print(df, ...)
  cat("\n")

  best <- which(apply(
    x$params, 1, function(z) identical(z, x$best_params)
  ))
  optimal_y <- x$stats[best, "mean"]

  cat(strwrap(paste0(
    "The best combination is: ",
    x$stats$combinations[best],
    " for a mean ", x$metric, " of ",
    format(optimal_y, digits = 3), ".\n"
  ), getOption("width")))
}
