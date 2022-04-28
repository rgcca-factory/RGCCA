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
#'   politic = Russett[, 6:11]
#' )
#' res <- rgcca_cv(blocks,
#'   response = 3, method = "rgcca",
#'   par_type = "sparsity",
#'   par_value = c(0.6, 0.75, 0.5),
#'   n_run = 2, n_cores = 1
#' )
#' plot(res)
#' \dontrun{
#' rgcca_cv(blocks,
#'   response = 3, par_type = "tau",
#'   par_value = c(0.6, 0.75, 0.5),
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
                     ...) {
  ### Try to retrieve parameters from a rgcca object
  if (!missing(blocks) & class(blocks) == "rgcca") {
    rgcca_res <- blocks
  }
  if (class(rgcca_res) == "rgcca") {
    message("All parameters were imported from a rgcca object.")
    scale_block <- rgcca_res$call$scale_block
    scale <- rgcca_res$call$scale
    scheme <- rgcca_res$call$scheme
    response <- rgcca_res$call$response
    tol <- rgcca_res$call$tol
    NA_method <- rgcca_res$call$NA_method
    bias <- rgcca_res$call$bias
    blocks <- rgcca_res$call$raw
    superblock <- rgcca_res$call$superblock
    tau <- rgcca_res$call$tau
    ncomp <- rgcca_res$call$ncomp
    sparsity <- rgcca_res$call$sparsity
  }

  ### Check parameters
  if (is.null(response)) {
    stop(paste0(
      "response is required for rgcca_cv (it is an integer ",
      "comprised between 1 and the number of blocks) "
    ))
  }
  if (validation == "loo") {
    k <- dim(blocks[[1]])[1]
    n_run <- 1
  }
  tmp <- check_prediction_model(prediction_model, blocks[[response]])

  check_integer("par_length", par_length)
  check_integer("n_run", n_run)
  check_integer("k", k, min = 2)
  match.arg(par_type, c("tau", "sparsity", "ncomp"))
  match.arg(validation, c("loo", "kfold"))

  ### Prepare parameters for line search
  if (method %in% c("sgcca", "spca", "spls") && (par_type == "tau")) {
    par_type <- "sparsity"
  } else if (par_type == "sparsity") method <- "sgcca"

  param <- set_parameter_grid(par_type, par_length, par_value, blocks, response)

  # Generate a warning if tau is null for a block that has more columns
  # than samples
  n <- NROW(blocks[[1]])
  overfitting_risk <- (param$par_type == "tau") &&
    any(vapply(seq_along(blocks), function(j) {
      NCOL(blocks[[j]]) > n && any(param$par_value[, j] == 0)
    }, FUN.VALUE = logical(1L)))
  if (overfitting_risk) {
    warning(
      "overfitting risk. Tau is zero for a block that has more columns ",
      "than rows, there is a high risk of overfitting for RGCCA."
    )
  }

  ### Start line search
  if (verbose) {
    message(paste("Cross-validation for", param$par_type, "in progress...\n"),
      appendLF = FALSE
    )
    pb <- txtProgressBar(max = nrow(param$par_value))
  }

  rgcca_args <- list(
    blocks = blocks, method = method, response = response,
    ncomp = ncomp, superblock = superblock,
    scale = scale, scale_block = scale_block,
    scheme = scheme, tol = tol, NA_method = NA_method,
    tau = tau, sparsity = sparsity, bias = bias,
    init = init
  )

  mat_cval <- lapply(seq(nrow(param$par_value)), function(i) {
    rgcca_args[[param$par_type]] <- param$par_value[i, ]
    rgcca_res <- do.call(rgcca, rgcca_args)

    res <- unlist(lapply(seq(n_run), function(n) {
      rgcca_cv_k(
        rgcca_res,
        validation = validation,
        prediction_model = tmp$prediction_model,
        k = k,
        n_cores = n_cores,
        complete = FALSE,
        verbose = FALSE,
        classification = tmp$classification,
        ...
      )$vec_scores
    }))

    if (verbose) setTxtProgressBar(pb, i)

    return(res)
  })
  mat_cval <- do.call(rbind, mat_cval)

  ### Format output
  cat("\n")

  connection <- matrix(0, nrow = length(blocks), ncol = length(blocks))
  connection[response, ] <- connection[, response] <- 1
  connection[response, response] <- 0

  call <- list(
    n_run = n_run,
    response = response,
    par_type = par_type,
    par_value = par_value,
    validation = validation,
    prediction_model = prediction_model,
    k = k,
    superblock = FALSE,
    scale = scale,
    scale_block = scale_block,
    tol = tol,
    scheme = scheme,
    NA_method = NA_method,
    blocks = blocks,
    tau = tau,
    sparsity = sparsity,
    ncomp = ncomp,
    connection = connection,
    init = init,
    bias = bias,
    method = method
  )

  rownames(param$par_value) <- seq(NROW(param$par_value))
  colnames(param$par_value) <- names(blocks)
  rownames(mat_cval) <- seq(NROW(mat_cval))

  best_param_idx <- which.min(apply(mat_cval, 1, mean))

  res <- list(
    cv = mat_cval,
    call = call,
    bestpenalties = param$par_value[best_param_idx, ],
    penalties = param$par_value
  )
  class(res) <- "cval"
  return(res)
}
