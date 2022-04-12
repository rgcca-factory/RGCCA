#' Cross-validation
#'
#' Uses cross-validation to evaluate predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams rgcca
#' @inheritParams bootstrap
#' @inheritParams plot_ind
#' @param classification A logical indicating if it is a classification task.
#' @noRd
rgcca_cv_k <- function(rgcca_res,
                       validation = "kfold",
                       prediction_model = "lm",
                       k = 5,
                       scale = NULL,
                       scale_block = NULL,
                       tol = 1e-8,
                       scheme = NULL,
                       NA_method = NULL,
                       method = NULL,
                       init = NULL,
                       bias = NULL,
                       connection = NULL,
                       ncomp = NULL,
                       tau = NULL,
                       sparsity = NULL,
                       n_cores = 1,
                       verbose = TRUE,
                       classification = FALSE,
                       ...) {
  ### Check parameters
  all_args <- names(environment())
  used_args <- c(
    names(match.call()), "validation", "prediction_model",
    "k", "n_cores", "verbose", "classification"
  )
  for (n in setdiff(all_args, used_args)) {
    assign(n, rgcca_res$call[[n]])
  }

  ### Compute cross validation
  blocks <- rgcca_res$call$raw

  if (validation == "loo") {
    v_inds <- seq(nrow(blocks[[1]]))
  } else {
    if (classification) {
      v_inds <- caret::createFolds(
        blocks[[rgcca_res$call$response]][, 1],
        k = k, list = TRUE,
        returnTrain = FALSE
      )
    } else {
      v_inds <- sample(nrow(blocks[[1]]))
      v_inds <- split(v_inds, seq(v_inds) %% k)
    }
  }

  res <- par_pblapply(
    v_inds, function(inds) {
      # Fit RGCCA on the training blocks
      res <- set_rgcca(rgcca_res,
        blocks = blocks,
        scale = scale,
        scale_block = scale_block,
        tol = tol,
        scheme = scheme,
        superblock = FALSE,
        inds = inds,
        NA_method = NA_method,
        response = rgcca_res$call$response,
        bias = bias,
        tau = tau,
        ncomp = ncomp,
        sparsity = sparsity,
      )

      # Remove columns of the validation blocks that had null variance in the
      # training blocks.
      column_sd_null <- remove_null_sd(
        lapply(blocks, function(x) x[-inds, , drop = FALSE])
      )$column_sd_null
      blocks_test <- lapply(seq_along(blocks), function(j) {
        if (length(column_sd_null[[j]]) > 0) {
          return(blocks[[j]][inds, -column_sd_null[[j]], drop = FALSE])
        }
        return(blocks[[j]][inds, , drop = FALSE])
      })
      names(blocks_test) <- names(blocks)

      # Evaluate RGCCA on the validation blocks
      res_pred <- rgcca_predict(
        res,
        blocks_test = blocks_test,
        prediction_model = prediction_model,
        response = rgcca_res$call$response,
        ...
      )
    },
    n_cores = n_cores, verbose = verbose
  )

  ### Structure outputs
  vec_scores <- vapply(res, function(x) x$score, FUN.VALUE = numeric(1))
  scores <- mean(unlist(lapply(res, function(x) x$score)), na.rm = TRUE)

  structure(
    list(
      scores = scores, vec_scores = vec_scores
    ),
    class = "cv"
  )
}
