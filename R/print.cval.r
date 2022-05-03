#' print.cval
#'
#' @inheritParams plot.cval
#' @param ... Further print options
#' @export
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca", par_type = "tau",
#'   par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1
#' )
#' print(res)
print.cval <- function(x, type = "quantile", ...) {
  summary_cval <- function(x, type = "quantile") {
    ymean <- apply(x$cv, 1, mean)
    ymed <- apply(x$cv, 1, median)

    switch(type,
      "quantile" = {
        middle_name <- "Median error"
        middle <- ymed
        lower <- apply(x$cv, 1, quantile, 0.05)
        upper <- apply(x$cv, 1, quantile, 0.95)
        low_lim <- "2.5%"
        up_lim <- "97.5%"
      },
      "sd" = {
        middle_name <- "Mean error"
        middle <- ymean
        lower <- middle - apply(x$cv, 1, sd)
        upper <- middle + apply(x$cv, 1, sd)
        low_lim <- "Mean - Sd"
        up_lim <- "Mean + Sd"
      },
      "stderr" = {
        middle_name <- "Mean error"
        middle <- ymean
        lower <- middle - apply(x$cv, 1, function(y) {
          sd(y) / sqrt(length(y))
        })
        upper <- middle + apply(x$cv, 1, function(y) {
          sd(y) / sqrt(length(y))
        })
        low_lim <- "Mean - Std Error"
        up_lim <- "Mean + Std Error"
      }
    )

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

    df <- data.frame(
      config = combinations, format(cbind(middle, lower, upper), digits = 3)
    )
    colnames(df) <- c("Tuning parameters", middle_name, low_lim, up_lim)
    return(df)
  }

  match.arg(type, c("sd", "stderr", "quantile"))

  ### Print parameters of the function
  print_call(x$call)

  penalties <- round(x$penalties, 3)
  rownames(penalties) <- seq(NROW(penalties))
  cat(fill = TRUE)
  cat(paste0("Tuning parameters (", x$call$par_type, ") used: "), fill = TRUE)
  print(penalties, quote = FALSE, ...)
  cat("\n")

  df <- summary_cval(x, type)

  cat(paste0(
    "Validation: ", x$call$validation,
    ifelse(x$call$validation == "kfold",
      paste0(
        " with ", x$call$k, " folds and ",
        x$call$n_run, " run(s))"
      ), ""
    )
  ), "\n")
  cat(paste("Prediction model:", x$call$prediction_model, "\n"))

  cat("\n")
  print(df)
  cat("\n")

  optimal_ind <- which(apply(
    x$penalties, 1, function(z) identical(z, x$bestpenalties)
  ))
  optimal_y <- df[optimal_ind, 2]

  cat(paste(
    "The best combination is:",
    paste(format(x$bestpenalties, digits = 3), collapse = " "),
    "for a mean CV error of",
    format(optimal_y, digits = 3)
  ), "\n", ...)
}
