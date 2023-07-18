#' @export
#' @rdname print
#' @order 2
print.rgcca_cv <- function(x, ...) {
  stopifnot(is(x, "rgcca_cv"))
  if (x$validation == "kfold") {
    cat(
      "RGCCA cross-validation model fitted on", nrow(x$params),
      "sets of parameters using", x$k, "folds."
    )
  } else {
    cat(
      "RGCCA cross-validation model fitted on", nrow(x$params),
      "sets of parameters using leave-one-out strategy."
    )
  }
  cat("\n")
  task <- ifelse(x$classification, "Classification", "Regression")
  cat(task, "was performed using a", x$prediction_model, "model.")
  cat("\n")
}
