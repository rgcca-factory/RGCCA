#' Tune RGCCA parameters by cross-validation
#'
#' This function is used to select automatically "sparsity", "tau" or "ncomp"
#' by cross-validation. This function only applies in a supervised setting,
#' and filling the response argument is therefore mandatory.
#'
#' @inheritParams rgcca
#' @inheritParams rgcca_predict
#' @inheritParams rgcca_bootstrap
#' @param par_type A character giving the parameter to tune among "sparsity",
#' "tau" or "ncomp".
#' @param par_value The parameter values to be tested, either NULL,
#' a numerical vector of size \eqn{J}{J}, or a matrix of size
#' par_length \eqn{\times J}{x J}.
#'
#' If par_value is NULL, up to par_length sets of parameters are generated
#' uniformly from
#' the minimum and maximum possible values of the parameter defined by par_type
#' for each block. Minimum possible values are 0 for tau,
#' \eqn{1/\textrm{sqrt}(p_j)}{1/sqrt(p_j)} for sparsity, and 1
#' for ncomp. Maximum possible values are 1 for tau and sparsity, and
#' \eqn{p_j}{p_j} for ncomp.
#'
#' If par_value is a vector, it overwrites the maximum values taken for the
#' range of generated parameters.
#'
#' If par_value is a matrix, par_value directly corresponds to the set of
#' tested parameters.
#' @param par_length An integer indicating the number of sets of candidate
#' parameters to be tested (if par_value is not a matrix).
#' @param k An integer giving the number of folds (if validation = 'kfold').
#' @param validation A string specifying the type of validation among "loo" and
#' "kfold". For small datasets (e.g. <30 samples), it is recommended to use a
#' loo (leave-one-out) procedure.
#' @param n_run An integer giving the number of Monte-Carlo Cross-Validation
#' (MCCV) to be run (if validation = 'kfold').
#' @export
#' @return A rgcca_cv object that can be printed and plotted.
#' @return  \item{k}{An integer giving the number of folds.}
#' @return  \item{n_run}{An integer giving the number of MCCV.}
#' @return  \item{opt}{A list containing some options of the
#' RGCCA model.}
#' @return  \item{metric}{A string indicating the metric used during the process
#' of cross-validation.}
#' @return \item{cv}{A matrix of dimension par_length x (k x n_run).
#' Each row of cv
#' corresponds to one set of candidate parameters. Each column of cv corresponds
#' to the cross-validated score of a specific fold in a specific run.}
#' @return \item{call}{A list of the input parameters of the RGCCA model.}
#' @return \item{par_type}{The type of parameter tuned (either "tau",
#' "sparsity", or "ncomp").}
#' @return \item{best_params}{The set of parameters that yields the best
#' cross-validated scores.}
#' @return \item{params}{A matrix reporting the sets of candidate parameters
#' used during the cross-validation process.}
#' @return \item{validation}{A string specifying the type of validation
#' (either "loo" or "kfold").}
#' @return \item{stats}{A data.frame containing various statistics (mean, sd,
#' median, first quartile, third quartile) of the cross-validated score for
#' each set of parameters that has been tested.}
#' @return \item{classification }{A boolean indicating if the model performs a
#' classification task.}
#' @return \item{prediction_model }{A string giving the model used for
#' prediction.}
#' @details
#' If the response block is univariate. The RGCCA components of each block
#' are used as input variables of the predictive model (specified by
#' "prediction_model") to predict the response block. The best combination of
#' parameters is the one with the best cross-validated score.
#' For multivariate response block, The RGCCA components of each block
#' are used as input variables of the predictive models (specified by
#' "prediction_model") to predict each column of the response block.
#' The cross-validated scores of each model are then averaged. The best
#' combination of parameters is the one with the best averaged cross-validated
#' score.
#' @examples
#'
#' # Cross_validation for classification
#'
#' set.seed(27) #favorite number
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, 1:3],
#'   industry = Russett[, 4:5],
#'   politic = as.factor(apply(Russett[, 9:11], 1, which.max))
#' )
#'
#' cv_out <- rgcca_cv(blocks, response = 3, method = "rgcca",
#'                    par_type = "tau",
#'                    par_length = 5,
#'                    prediction_model = "lda", #caret::modelLookup()
#'                    metric = "Accuracy",
#'                    k=3, n_run = 3,
#'                    verbose = TRUE)
#'
#'
#' print(cv_out)
#' plot(cv_out)
#'
#' # A fitted cval object is given as output of the rgcca() function
#'
#' fit_opt = rgcca(cv_out)
#' \dontrun{
#' # Cross_validation for regression
#'
#' set.seed(27) #favorite number
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, 1:3],
#'   industry = Russett[, 4:5],
#'   politic =  Russett[, 6:8]
#' )
#'
#' cv_out <- rgcca_cv(blocks, response = 3, method = "rgcca",
#'                    par_type = "tau",
#'                    par_value = c(0.6, 0.75, 0.8),
#'                    prediction_model = "lm", #caret::modelLookup()
#'                    metric = "RMSE",
#'                    k=3, n_run = 5,
#'                    verbose = TRUE)
#'
#' print(cv_out)
#' plot(cv_out)
#'
#' fit_opt = rgcca(cv_out)
#'
#'
#'  data("ge_cgh_locIGR", package = "gliomaData")
#'  blocks <- ge_cgh_locIGR$multiblocks
#'  Loc <- factor(ge_cgh_locIGR$y)
#'  levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#'  blocks[[3]] <- Loc
#'  set.seed(27) # favorite number
#'
#'   cv_out = rgcca_cv(blocks, response = 3,
#'                    ncomp = 1,
#'                    prediction_model = "glmnet",
#'                    family = "multinomial", lambda = .001,
#'                    par_type = "sparsity",
#'                    par_value = c(.071, .2, 1),
#'                    metric = "Balanced_Accuracy",
#'                    n_cores = 2,
#'  )
#'
#'  print(cv_out)
#'  plot(cv_out, display_order = FALSE)
#'
#'   cv_out = rgcca_cv(blocks, response = 3,
#'                    ncomp = 1,
#'                    prediction_model = "glmnet",
#'                    family = "multinomial", lambda = .001,
#'                    par_type = "ncomp",
#'                    par_value = c(5, 5, 1),
#'                    metric = "Balanced_Accuracy",
#'                    n_cores = 2,
#'  )
#'
#'  print(cv_out)
#'  plot(cv_out, display_order = FALSE)
#' }
rgcca_cv <- function(blocks,
                     connection = NULL,
                     method = "rgcca",
                     response = NULL,
                     par_type = "tau",
                     par_value = NULL,
                     par_length = 10,
                     validation = "kfold",
                     prediction_model = "lm",
                     metric = NULL,
                     k = 5,
                     n_run = 1,
                     n_cores = 1,
                     quiet = TRUE,
                     superblock = FALSE,
                     scale = TRUE,
                     scale_block = TRUE,
                     tol = 1e-8,
                     scheme = "factorial",
                     NA_method = "na.ignore",
                     rgcca_res = NULL,
                     tau = 1,
                     ncomp = 1,
                     sparsity = 1,
                     init = "svd",
                     bias = TRUE,
                     verbose = TRUE,
                     n_iter_max = 1000,
                     comp_orth = TRUE,
                     rank = 1,
                     mode_orth = 1,
                     separable = TRUE,
                     ...) {
  ### Try to retrieve parameters from a rgcca object
  rgcca_args <- as.list(environment())
  tmp <- get_rgcca_args(blocks, rgcca_args)
  opt <- tmp$opt
  rgcca_args <- tmp$rgcca_args

  if (is.null(rgcca_args$response)) {
    stop(paste0(
      "response is required for rgcca_cv (it is an integer ",
      "comprised between 1 and the number of blocks) "
    ))
  }

  ### Check parameters
  if (validation == "loo") {
    k <- NROW(rgcca_args$blocks[[1]])
    n_run <- 1
  }
  model <- check_prediction_model(
    prediction_model, rgcca_args$blocks[[rgcca_args$response]],
    missing(prediction_model)
  )

  check_integer("par_length", par_length)
  check_integer("n_run", n_run)
  check_integer("k", k, min = 2)
  match.arg(par_type, c("tau", "sparsity", "ncomp"))
  match.arg(validation, c("loo", "kfold"))

  default_metric <- ifelse(model$classification, "Accuracy", "RMSE")
  metric <- ifelse(is.null(metric), default_metric, metric)
  available_metrics <- get_available_metrics(model$classification)
  metric <- match.arg(metric, available_metrics)

  ### Prepare parameters for line search
  if (
    rgcca_args$method %in% sparse_methods() && (par_type == "tau")
  ) {
    par_type <- "sparsity"
  } else if (par_type == "sparsity") {
    rgcca_args$method <- "sgcca"
    opt$param <- "sparsity"
  }

  param <- set_parameter_grid(
    par_type, par_length, par_value, rgcca_args$blocks,
    rgcca_args[[par_type]], rgcca_args$response, FALSE, opt$disjunction
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

  ### Create folds
  idx <- seq_len(NROW(rgcca_args$blocks[[1]]))
  if (validation == "loo") {
    v_inds <- idx
  } else {
    if (model$classification) {
      v_inds <- Reduce("c", lapply(seq_len(n_run), function(j) {
        folds <- caret::createFolds(
          rgcca_args$blocks[[rgcca_args$response]][, 1],
          k = k, list = TRUE,
          returnTrain = FALSE
        )
        # If there are NA in the response block, caret creates an extra fold
        # with all the NA elements. We split them across the other folds
        if (length(folds) > k) {
          extra_folds <- split(folds[[1]], seq_along(folds[[1]]) %% k)
          folds <- folds[-1]
          folds <- Map(append, folds, extra_folds)
        }
        return(folds)
      }))
    } else {
      v_inds <- Reduce("c", lapply(seq_len(n_run), function(j) {
        x <- sample(idx)
        split(x, seq(x) %% k)
      }))
    }
  }

  ### Compute cross validation
  idx <- seq_len(NROW(param$par_value) * length(v_inds))
  W <- par_pblapply(idx, function(n) {
    i <- (n - 1) %/% length(v_inds) + 1
    j <- (n - 1) %% length(v_inds) + 1
    rgcca_cv_k(
      rgcca_args,
      inds = v_inds[[j]],
      metric = metric,
      par_type = param$par_type,
      par_value = param$par_value[i, ],
      prediction_model = model$prediction_model,
      ...
    )
  }, n_cores = n_cores, verbose = verbose)

  W <- matrix(unlist(W), nrow = NROW(param$par_value), byrow = TRUE)

  if (model$classification) {
    W[is.na(W)] <- 0
  }

  ### Format output
  rownames(param$par_value) <- seq_len(NROW(param$par_value))
  colnames(param$par_value) <- names(rgcca_args$blocks)
  rownames(W) <- seq_len(NROW(W))

  best_param_idx <- ifelse(
    model$classification,
    which.max(apply(W, 1, mean, na.rm = TRUE)),
    which.min(apply(W, 1, mean, na.rm = TRUE))
  )

  # Compute statistics
  combinations <- format_combinations(param$par_value)

  stats <- data.frame(
    combinations,
    mean = apply(W, 1, mean),
    sd = apply(W, 1, sd),
    median = apply(W, 1, median),
    Q1 = apply(W, 1, quantile, 0.25),
    Q3 = apply(W, 1, quantile, 0.75)
  )

  res <- list(
    k = k,
    cv = W,
    opt = opt,
    call = rgcca_args,
    stats = stats,
    n_run = n_run,
    metric = metric,
    par_type = param$par_type,
    params = param$par_value,
    validation = validation,
    best_params = param$par_value[best_param_idx, ],
    classification = model$classification,
    prediction_model = model$model_name
  )
  class(res) <- "rgcca_cv"
  return(res)
}
