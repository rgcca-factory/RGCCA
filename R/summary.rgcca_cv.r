#' @export
#' @rdname summary
#' @order 2
summary.rgcca_cv <- function(object, type = c("sd", "quantile"), ...) {
  stopifnot(is(object, "rgcca_cv"))
  type <- type[1]
  type <- match.arg(type, c("sd", "quantile"))

  ### Print parameters of the function
  
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
  optimal_y <- object$stats[object$best_param_idx, "mean"]


  cat(strwrap(paste0(
    "The best combination is: ",
   list(round(unlist(object$best_params), 3)),
    " for a mean ", object$metric, " of ",
    format(optimal_y, digits = 3)
  ), getOption("width")))
  cat("\n")
}
