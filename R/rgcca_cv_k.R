#' Cross-validation
#'
#' Uses cross-validation to evaluate predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams rgcca
#' @inheritParams bootstrap
#' @inheritParams plot_ind
#' @param k An integer giving the number of folds (if validation = 'kfold').
#' @param validation A character for the type of validation among "loo",
#' "kfold", "test".
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' rgcca_out <- rgcca(blocks, response = 3, superblock = FALSE)
#' res <- rgcca_cv_k(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
#' rgcca_cv_k(rgcca_out, n_cores = 1)
#' @export
#' @seealso \link{rgcca}, \link{rgcca_predict}, \link{plot.predict}
rgcca_cv_k <- function(rgcca_res,
                       validation = "kfold",
                       task = "regression",
                       prediction_model = "lm",
                       X_scaled = TRUE,
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
                       n_cores = parallel::detectCores() - 1,
                       verbose = TRUE,
                       ...) {
  ### Check parameters
  stopifnot(is(rgcca_res, "rgcca"))
  if (is.null(rgcca_res$call$response)) {
    stop_rgcca(
      "missing response block. A model with a response block must be ",
      "used to apply rgcca_cv_k."
    )
  }
  match.arg(validation, c("loo", "test", "kfold"))

  all_args <- names(environment())
  used_args <- c(
    names(match.call()), "validation", "task", "prediction_model",
    "X_scaled", "k", "n_cores", "verbose"
  )
  for (n in setdiff(all_args, used_args)) {
    assign(n, rgcca_res$call[[n]])
  }

  check_integer("k", k, min = 2)
  check_integer("n_cores", n_cores, min = 0)

  ### Compute cross validation
  blocks <- rgcca_res$call$raw

  if (validation == "loo") {
    v_inds <- seq(nrow(blocks[[1]]))
  } else {
    v_inds <- sample(nrow(blocks[[1]]))
    v_inds <- split(v_inds, seq(v_inds) %% k)
  }

  scores <- par_pblapply(
    v_inds, function(inds) {
      res <- set_rgcca(rgcca_res,
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
      # TODO: decide what to do with null variance variables, is it ok to just
      # keep them (by passing init = FALSE to check_blocks?)
      res_pred <- rgcca_predict(
        res,
        X = lapply(blocks, function(x) x[inds, , drop = FALSE]),
        task = task,
        prediction_model = prediction_model,
        block_to_predict = rgcca_res$call$response,
        X_scaled = FALSE
      )
    },
    n_cores = n_cores, verbose = verbose
  )

  list_rgcca <- lapply(scores, function(x) {
    return(x$rgcca_res)
  })
  list_pred <- lapply(scores, function(x) {
    return(x$pred)
  })
  list_scores <- sapply(scores, function(x) x$score)
  list_res <- lapply(scores, function(x) {
    return(x$res)
  })
  list_class_fit <- lapply(scores, function(x) {
    return(x$class.fit)
  })

  # concatenation of each test set to provide predictions for each block
  preds <- lapply(
    seq_along(rgcca_res$call$blocks),
    function(x) {
      Reduce(
        function(y, z) {
          rbind(y, z$pred[[x]])
        },
        scores,
        init = NULL
      )
    }
  )

  names(preds) <- names(rgcca_res$call$blocks)

  for (x in seq_along(preds)) {
    preds[[x]] <- preds[[x]][rownames(blocks[[x]]), , drop = FALSE]
  }

  scores <- mean(unlist(lapply(scores, function(x) x$score)), na.rm = T)

  structure(
    list(
      scores = scores, preds = preds,
      rgcca_res = rgcca_res,
      list_scores = list_scores,
      list_pred = list_pred, list_rgcca = list_rgcca,
      list_class = list_class_fit, list_res = list_res
    ),
    class = "cv"
  )
}
