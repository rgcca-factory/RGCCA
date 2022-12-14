#' Cross-validation
#'
#' Uses cross-validation to evaluate predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams rgcca
#' @inheritParams bootstrap
#' @noRd
rgcca_cv_k <- function(rgcca_args, inds, prediction_model,
                       par_type, par_value, metric, ...) {
  rgcca_args[[par_type]] <- par_value

  blocks <- rgcca_args[["blocks"]]

  rgcca_args[["blocks"]] <- lapply(
    blocks, function(x) x[-inds, , drop = FALSE]
  )
  # Fit RGCCA on the training blocks
  res <- do.call(rgcca, rgcca_args)

  # Evaluate RGCCA on the validation blocks
  blocks_test <- lapply(seq_along(blocks), function(j) {
    x <- blocks[[j]][inds, , drop = FALSE]
    colnames(x) <- colnames(res$call$blocks[[j]])
    return(x)
  })
  names(blocks_test) <- names(res$blocks)

  return(rgcca_predict(
    res,
    metric = metric,
    blocks_test = blocks_test,
    prediction_model = prediction_model,
    ...
  )$score)
}
