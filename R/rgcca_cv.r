#' Tune RGCCA parameters in 'supervised' mode with cross-validation
#'
#' Tune the sparsity coefficient (if the model is sparse) or tau
#' (otherwise) in a supervised approach by estimating by crossvalidation the
#' predictive quality of the models.
#' In this purpose, the samples are divided into k folds where the model will
#' be tested on each fold and trained on the others. For small datasets
#' (<30 samples), it is recommended to use as many folds as there are
#' individuals (leave-one-out; loo).
#' @inheritParams rgcca
#' @inheritParams rgcca_predict
#' @inheritParams bootstrap
#' @param par_type A character giving the parameter to tune among "sparsity"
#' or "tau".
#' @param par_value A matrix (n*p, with p the number of blocks and n the number
#' of combinations to be tested), a vector (of p length) or a numeric value
#' giving sets of penalties (tau for RGCCA, sparsity for SGCCA) to be tested,
#' one row by combination. By default, it takes 10 sets between min values (0
#'  for RGCCA and $1/sqrt(ncol)$ for SGCCA) and 1.
#' @param par_length An integer indicating the number of sets of parameters to
#' be tested (if par_value = NULL). The parameters are uniformly distributed.
#' @param k An integer giving the number of folds (if validation = 'kfold').
#' @param validation A character for the type of validation among "loo",
#' "kfold".
#' @param n_run An integer giving the number of cross-validations to be run
#' (if validation = 'kfold').
#' @export
#' @return \item{cv}{A matrix giving the root-mean-square error (RMSE) between
#' the predicted R/SGCCA and the observed R/SGCCA for each combination and each
#' prediction (n_prediction = n_samples for validation = 'loo';
#' n_prediction = 'k' * 'n_run' for validation = 'kfold').}
#' @return \item{call}{A list of the input parameters}
#' @return \item{bestpenalties}{Penalties giving the best RMSE for each blocks
#' (for regression) or the best proportion of wrong predictions
#' (for classification)}
#' @return \item{penalties}{A matrix giving, for each blocks, the penalty
#' combinations (tau or sparsity)}
#' @details
#' At each round of cross-validation, for each
#' variable, a predictive model of the first RGCCA component of each block
#' (calculated on the training set) is constructed.
#' Then the Root Mean Square of Errors (RMSE) or the Accuracy of the model
#' is computed on the testing dataset. Finally, the metrics are averaged on the
#' different folds.
#' The best combination of parameters is the one where the average of RMSE on
#' the testing datasets is the lowest or the accuracy is the highest.
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:8]
#' )
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca",
#'   par_type = "sparsity",
#'   par_value = c(0.6, 0.75, 0.8),
#'   n_run = 2, n_cores = 1
#' )
#' plot(res)
#' \dontrun{
#' rgcca_cv(blocks,
#'   response = 3, par_type = "tau",
#'   par_value = c(0.6, 0.75, 0.8),
#'   n_run = 2, n_cores = 1
#' )$bestpenalties
#'
#' rgcca_cv(blocks,
#'   response = 3, par_type = "sparsity",
#'   par_value = 0.8, n_run = 2, n_cores = 1
#' )
#'
#' rgcca_cv(blocks,
#'   response = 3, par_type = "tau",
#'   par_value = 0.8, n_run = 2, n_cores = 1
#' )
#' }
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
rgcca_cv <- function(blocks,
                     method = "rgcca",
                     response = NULL,
                     par_type = "tau",
                     par_value = NULL,
                     par_length = 10,
                     validation = "kfold",
                     prediction_model = "lm",
                     k = 5,
                     n_run = 1,
                     n_cores = 1,
                     quiet = TRUE,
                     superblock = FALSE,
                     scale = TRUE,
                     scale_block = TRUE,
                     tol = 1e-8,
                     scheme = "factorial",
                     NA_method = "nipals",
                     rgcca_res = NULL,
                     tau = rep(1, length(blocks)),
                     ncomp = rep(1, length(blocks)),
                     sparsity = rep(1, length(blocks)),
                     init = "svd",
                     bias = TRUE,
                     verbose = TRUE,
                     n_iter_max = 1000,
                     score = NULL,
                     ...) {
  ### Try to retrieve parameters from a rgcca object
  rgcca_args <- as.list(environment())
  rgcca_args <- get_rgcca_args(blocks, rgcca_args)
  rgcca_args$quiet <- TRUE
  rgcca_args$verbose <- FALSE

  rgcca_args$blocks <- check_blocks(
    rgcca_args$blocks,
    add_NAlines = TRUE, quiet = rgcca_args$quiet,
    response = rgcca_args$response
  )

  ### Check parameters
  if (is.null(rgcca_args$response)) {
    stop(paste0(
      "response is required for rgcca_cv (it is an integer ",
      "comprised between 1 and the number of blocks) "
    ))
  }
  check_blockx("response", rgcca_args$response, rgcca_args$blocks)
  if (validation == "loo") {
    k <- dim(rgcca_args$blocks[[1]])[1]
    n_run <- 1
  }
  model <- check_prediction_model(
    prediction_model, rgcca_args$blocks[[rgcca_args$response]]
  )

  check_integer("par_length", par_length)
  check_integer("n_run", n_run)
  check_integer("k", k, min = 2)
  match.arg(par_type, c("tau", "sparsity", "ncomp"))
  match.arg(validation, c("loo", "kfold"))

  default_score <- ifelse(model$classification, "Accuracy", "RMSE")
  score <- ifelse(is.null(score), default_score, score)
  available_scores <- get_available_scores(model$classification)
  score <- match.arg(score, available_scores)

  ### Set connection matrix
  connection <- matrix(
    0, nrow = length(rgcca_args$blocks), ncol = length(rgcca_args$blocks)
  )
  connection[rgcca_args$response, ] <- connection[, rgcca_args$response] <- 1
  connection[rgcca_args$response, rgcca_args$response] <- 0
  rgcca_args$connection <- connection

  ### Prepare parameters for line search
  if (
    rgcca_args$method %in% c("sgcca", "spca", "spls") && (par_type == "tau")
  ) {
    par_type <- "sparsity"
  } else if (par_type == "sparsity") {
    rgcca_args$method <- "sgcca"
  }

  param <- set_parameter_grid(
    par_type, par_length, par_value, rgcca_args$blocks, rgcca_args$response
  )

  # Generate a warning if tau has not been fully specified for a block that
  # has more columns than samples and remove tau = 0 configuration
  n <- NROW(rgcca_args$blocks[[1]])
  overfitting_risk <- (param$par_type == "tau") && is.null(dim(par_value)) &&
    any(vapply(
      seq_along(rgcca_args$blocks),
      function(j) NCOL(rgcca_args$blocks[[j]]) > n,
      FUN.VALUE = logical(1L)
    ))
  if (overfitting_risk) {
    param$par_value <- param$par_value[-nrow(param$par_value), ]
    warning(
      "overfitting risk. A block has more columns than rows, so the ",
      "configuration with tau = 0 has been removed."
    )
  }

  ### Start line search
  if (validation == "loo") {
    v_inds <- seq_len(NROW(rgcca_args$blocks[[1]]))
  } else {
    if (model$classification) {
      v_inds <- Reduce("c", lapply(seq_len(n_run), function(j) {
        caret::createFolds(
          rgcca_args$blocks[[rgcca_args$response]][, 1],
          k = k, list = TRUE,
          returnTrain = FALSE
        )
      }))
    } else {
      v_inds <- Reduce("c", lapply(seq_len(n_run), function(j) {
        x <- sample(nrow(rgcca_args$blocks[[1]]))
        split(x, seq(x) %% k)
      }))
    }
  }

  idx <- seq_len(NROW(param$par_value) * length(v_inds))
  W <- par_pblapply(idx, function(n) {
    i <- (n - 1) %/% length(v_inds) + 1
    j <- (n - 1) %% length(v_inds) + 1
    rgcca_cv_k(
      rgcca_args,
      inds = v_inds[[j]],
      score = score,
      par_type = param$par_type,
      par_value = param$par_value[i, ],
      prediction_model = model$prediction_model,
      ...
    )
  }, n_cores = n_cores, verbose = verbose)

  param$par_value[, rgcca_args$response] <- do.call(
    rbind,
    lapply(W[seq(1, length(idx), by = length(v_inds))], "[[", "par_value")
  )
  W <- unlist(lapply(W, "[[", "score"))
  W <- matrix(W, nrow = NROW(param$par_value), byrow = TRUE)

  ### Format output
  rownames(param$par_value) <- seq_len(NROW(param$par_value))
  colnames(param$par_value) <- names(rgcca_args$blocks)
  rownames(W) <- seq_len(NROW(W))

  best_param_idx <- ifelse(
    model$classification,
    which.max(apply(W, 1, mean)),
    which.min(apply(W, 1, mean))
  )

  res <- list(
    k = k,
    cv = W,
    call = rgcca_args,
    n_run = n_run,
    score = score,
    par_type = param$par_type,
    penalties = param$par_value,
    validation = validation,
    bestpenalties = param$par_value[best_param_idx, ],
    prediction_model = prediction_model
  )
  class(res) <- "cval"
  return(res)
}
