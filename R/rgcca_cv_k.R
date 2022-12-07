#' Cross-validation
#'
#' Uses cross-validation to evaluate predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams rgcca
#' @inheritParams bootstrap
#' @noRd
rgcca_cv_k <- function(rgcca_args, inds, prediction_model,
                       par_type, par_value, ...) {
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
  score <- rgcca_predict(
    res,
    blocks_test = blocks_test,
    prediction_model = prediction_model,
    response = rgcca_args$response,
    ...
  )$score

  ### Structure outputs
  structure(
    list(
      score = score, par_value = res$call[[par_type]][res$call$response]
    ),
    class = "cv"
  )
}
