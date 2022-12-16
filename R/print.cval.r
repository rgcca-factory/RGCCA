#' print.cval
#'
#' @inheritParams plot.cval
#' @param ... Further print options
#' @return none
#' @export
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:8]
#' )
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca", par_type = "tau",
#'   par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1
#' )
#' print(res)
print.cval <- function(x, type = "sd", ...) {
  type <- match.arg(type, c("sd", "quantile"))

  ### Print parameters of the function
  print_call(x$call)

  penalties <- round(x$penalties, 3)
  rownames(penalties) <- seq_len(NROW(penalties))
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", x$par_type, ") used: "), fill = TRUE)
  print(penalties, quote = FALSE, ...)
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
  print(df)
  cat("\n")

  optimal_ind <- which(apply(
    x$penalties, 1, function(z) identical(z, x$bestpenalties)
  ))
  optimal_y <- x$stats[optimal_ind, "mean"]

  cat(paste(
    "The best combination is:",
    paste(format(x$bestpenalties, digits = 3), collapse = " "),
    "for a mean", x$metric, "of",
    format(optimal_y, digits = 3)
  ), "\n", ...)
}
