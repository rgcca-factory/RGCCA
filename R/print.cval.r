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
print.cval <- function(x, bars = "quantile", ...) {
  summary_cval <- function(x, bars = "quantile") {
    mat_cval <- x$cv
    mean_b <- apply(mat_cval, 1, mean)

    match.arg(bars, c("sd", "stderr", "quantile"))
    if (bars != "none" && dim(mat_cval)[2] < 3) {
      bars == "none"
      warning(paste0(
        "Standard deviations can not be calculated with less ",
        "than 3 columns in x"
      ))
    }
    if (bars != "none") {
      if (bars == "quantile") {
        inf_b <- apply(mat_cval, 1, function(y) {
          return(quantile(y, 0.025))
        })
        sup_b <- apply(mat_cval, 1, function(y) {
          return(quantile(y, 0.975))
        })
        lowlim <- "2.5%"
        uplim <- "97.5%"
      }
      if (bars == "sd") {
        inf_b <- mean_b - apply(mat_cval, 1, sd)
        sup_b <- mean_b + apply(mat_cval, 1, sd)
        lowlim <- "Mean - Sd"
        uplim <- "Mean + Sd"
      }
      if (bars == "stderr") {
        inf_b <- mean_b - apply(mat_cval, 1, function(y) {
          sd(y) / sqrt(length(y))
        })
        sup_b <- mean_b + apply(mat_cval, 1, function(y) {
          sd(y) / sqrt(length(y))
        })
        lowlim <- "Mean - Std Error"
        uplim <- "Mean + Std Error"
      }
    }
    df <- round(data.frame(
      config = seq(nrow(mat_cval)), mean = mean_b,
      inf = inf_b, sup = sup_b
    ), 3)
    if (x$call$type_cv == "regression") {
      colnames(df) <- c("Combination", "Mean RMSE", lowlim, uplim)
    }
    if (x$call$type_cv == "classification") {
      colnames(df) <- c(
        "Combination", "Mean Error Prediction Rate",
        lowlim, uplim
      )
    }
    return(df)
  }


  cat("Call: ")
  names_call <- c(
    "type_cv", "n_run", "NA_method", "tol", "scale",
    "scale_block"
  )
  char_to_print <- ""
  for (name in names_call) {
    if (name == "ncomp") {
      if (length(x$call$ncomp) > 1) {
        value <- (paste(x$call$ncomp, sep = "", collapse = ","))
        value <- paste0("c(", value, ")")
      }
    }
    if (name != "ncomp") {
      value <- x$call[[name]]
    }
    quo <- ifelse(is.character(value) & name != "ncomp", "'", "")
    vir <- ifelse(name == names_call[length(names_call)], "", ", ")
    char_to_print <- paste(char_to_print, name, "=", quo, value, quo, vir,
      collapse = "", sep = ""
    )
  }
  cat(char_to_print)
  cat("\n")

  c1s <- round(x$penalties, 4)
  rownames(c1s) <- seq(NROW(c1s))
  cat(fill = TRUE)
  cat("Tuning parameters used: ", fill = TRUE)
  print(c1s, quote = FALSE, ...)
  cat("\n")

  df <- summary_cval(x, bars)
  colname_for_optimal <- ifelse(x$call$type_cv == "regression", "Mean RMSE",
    "Mean Error Prediction Rate"
  )
  optimal_ind <- which.min(df[, colname_for_optimal])
  optimal_y <- df[optimal_ind, colname_for_optimal]
  cat(paste0(nrow(x$cv), " configurations were tested. \n"))

  cat(paste0(
    "Validation: ", x$call$validation,
    ifelse(x$call$validation == "kfold",
      paste0(
        " with ", x$call$k, " folds and ",
        x$call$n_run, " run(s))"
      )
    )
  ), "\n")

  cat("\n")
  print(df)
  cat("\n")
  if (x$call$type_cv == "regression") {
    cat(paste(
      "The best combination was:",
      paste(round(x$bestpenalties, digits = 3), collapse = " "),
      "for a mean CV criterion (RMSE) of ",
      round(optimal_y, digits = 2)
    ), "\n", ...)
  }
  if (x$call$type_cv == "classification") {
    cat(paste(
      "The best combination was:",
      paste(round(x$bestpenalties, digits = 3), collapse = " "),
      "for a mean rate of false predictions of ",
      round(optimal_y, digits = 2)
    ), "\n", ...)
  }
}
