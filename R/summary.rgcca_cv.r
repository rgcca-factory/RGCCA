#' @export
#' @rdname summary
#' @order 2
summary.rgcca_cv <- function(object, type = c("sd", "quantile"), ...) {
  stopifnot(is(object, "rgcca_cv"))
  type <- type[1]
  type <- match.arg(type, c("sd", "quantile"))

  ### Print parameters of the function
  print_call(object$call)

  params <- round(object$params, 3)
  rownames(params) <- seq_len(NROW(params))
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", object$par_type, ") used: "), fill = TRUE)
  print(params, quote = FALSE, ...)
  cat("\n")

  cat(paste0(
    "Validation: ", object$validation,
    ifelse(object$validation == "kfold",
      paste0(
        " with ", object$k, " folds and ",
        object$n_run, " run(s))"
      ), ""
    )
  ), "\n")
  cat(paste("Prediction model:", object$prediction_model, "\n"))
  cat("\n")

  df <- format(object$stats, digits = 3)
  if (type == "quantile") {
    df <- df[, c(1, 4, 5, 6)]
    colnames(df) <- c(
      "Tuning parameters", paste("Median", object$metric), "Q1", "Q3"
    )
  } else {
    df <- df[, c(1, 2, 3)]
    colnames(df) <- c("Tuning parameters", paste("Mean", object$metric), "Sd")
  }
  print(df, ...)
  cat("\n")

  best <- which(apply(
    object$params, 1, function(z) identical(z, object$best_params)
  ))
  optimal_y <- object$stats[best, "mean"]

  cat(strwrap(paste0(
    "The best combination is: ",
    object$stats$combinations[best],
    " for a mean ", object$metric, " of ",
    format(optimal_y, digits = 3)
  ), getOption("width")))
  cat("\n")
}
