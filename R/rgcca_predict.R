#' Make predictions using RGCCA
#'
#' This function aims to make predictions combining a fitted RGCCA object
#' and a prediction model for classification or regression.
#'
#' @inheritParams rgcca_transform
#' @param blocks_test A list of test blocks from which we aim to predict the
#' associated response block. If the test response block is present among
#' blocks_test, metrics are computed by comparing the predictions and the
#' true values.
#' @param prediction_model A string giving the model used for prediction.
#' Please see caret::modelLookup() for a list of the available models.
#' @param metric A string indicating the metric of interest.
#' It should be one of the following scores:
#'
#' For classification: "Accuracy", "Kappa", "F1", "Sensitivity", "Specificity",
#' "Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", "Detection_Rate",
#' "Balanced_Accuracy".
#'
#' For regression: "RMSE", "MAE".
#' @param ... Additional parameters to be passed to prediction_model.
#' @return A list containing the following elements:
#' @return \item{score}{The score obtained on the testing block. NA if the test
#' block is missing.}
#' @return \item{model}{A list of the models trained using caret to make the
#' predictions and compute the scores.}
#' @return \item{probs}{A list of data.frames with the class probabilities
#' of the test and train response blocks predicted by the prediction
#' model. If the prediction model does not compute class probabilities, the
#' data.frames are empty.}
#' @return \item{metric}{A list of data.frames containing the scores obtained
#' on the training and testing sets.}
#' @return \item{confusion}{A list containing NA for regression tasks.
#' Otherwise, the confusion summary produced by caret for train and test.}
#' @return \item{projection}{A list of matrices containing the projections
#' of the test blocks using the canonical components from the fitted RGCCA
#' object. The response block is not projected.}
#' @return \item{prediction}{A list of data.frames with the predictions
#' of the test and train response blocks.}
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, 1:3],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:8]
#' )
#' X_train <- lapply(blocks, function(x) x[seq(1, 30), ])
#' X_test <- lapply(blocks, function(x) x[seq(31, 47), ])
#' fit <- rgcca(X_train,
#'   tau = 1, ncomp = c(3, 2, 3), response = 3
#' )
#' res <- rgcca_predict(fit, X_test)
#' @importFrom caret train trainControl confusionMatrix
#' @importFrom caret multiClassSummary postResample
#' @importFrom stats predict
#' @export
rgcca_predict <- function(rgcca_res,
                          blocks_test = rgcca_res$call$blocks,
                          prediction_model = "lm",
                          metric = NULL,
                          ...) {
  ### Check input parameters
  if (is.null(names(blocks_test))) {
    stop_rgcca("Please provide names for blocks_test.")
  }

  response <- rgcca_res$call$response
  if (is.null(response)) {
    stop_rgcca("RGCCA must use a response block.")
  }
  test_idx <- match(names(rgcca_res$blocks)[response], names(blocks_test))

  ### Check that response is among both training and test blocks
  if (is.na(test_idx)) {
    no_y_test <- TRUE
    n_test <- NROW(blocks_test[[1]])
    test_idx <- length(blocks_test) + 1
    blocks_test[[names(rgcca_res$blocks)[response]]] <- matrix(
      rnorm(n_test),
      nrow = n_test,
      ncol = NCOL(rgcca_res$call$blocks[[response]])
    )
  } else {
    no_y_test <- FALSE
  }

  tmp <- check_prediction_model(
    prediction_model, rgcca_res$call$blocks[[response]],
    missing(prediction_model)
  )
  prediction_model <- tmp$prediction_model
  classification <- tmp$classification

  default_metric <- ifelse(classification, "Accuracy", "RMSE")
  metric <- ifelse(is.null(metric), default_metric, metric)
  available_metrics <- get_available_metrics(classification)
  metric <- match.arg(metric, available_metrics)

  ### Get train and test target (if present)
  y_train <- rgcca_res$call$blocks[[response]]
  y_test <- as.data.frame(blocks_test[[test_idx]])

  if (any(dim(y_test)[-1] != dim(y_train)[-1])) {
    stop_rgcca(
      "Dimensions of response do not match between",
      " train and test blocks."
    )
  }

  ### Get projected train and test data
  projection <- rgcca_transform(rgcca_res, blocks_test[-test_idx])
  X_train <- rgcca_res$Y[names(projection)]
  X_train <- reformat_projection(X_train)
  X_test <- reformat_projection(projection)

  # Keep same lines in X_train and y_train
  y_train <- as.data.frame(subset_rows(y_train, rownames(X_train)))

  # Test that in classification, variables are not constant within groups
  if (classification) {
    groups <- split(X_train, y_train[, 1])
    is_constant <- unlist(lapply(groups, function(g) {
      apply(g, 2, function(x) {
        isTRUE(all.equal(
          x, rep(x[1], length(x)),
          check.attributes = FALSE
        ))
      })
    }))
    if (any(is_constant)) {
      warning(
        "overfitting risk. The RGCCA method led to projected blocks ",
        "that are constant within groups. Try to regularize the model",
        " by increasing tau or decreasing sparsity."
      )
    }
  }

  ### Train prediction model and predict results on X_test
  results <- lapply(
    seq_len(NCOL(y_train)), function(j) {
      core_prediction(
        prediction_model, X_train, X_test,
        y_train[, j], y_test[, j], metric,
        classification, no_y_test, ...
      )
    }
  )
  names(results) <- colnames(y_train)

  prediction <- lapply(c("train", "test"), function(mode) {
    as.data.frame(lapply(results, function(res) {
      res[["prediction"]][[mode]][, "pred"]
    }))
  })

  metric <- lapply(c("train", "test"), function(mode) {
    as.data.frame(lapply(results, function(res) {
      res[["metric"]][[mode]]
    }))
  })

  probs <- lapply(c("train", "test"), function(mode) {
    as.data.frame(lapply(results, function(res) {
      res[["probs"]][[mode]]
    }))
  })

  confusion <- results[[1]]$confusion

  names(prediction) <- names(metric) <- names(probs) <- c("train", "test")

  model <- lapply(results, "[[", "model")
  score <- mean(unlist(lapply(results, "[[", "score")), na.rm = TRUE)
  names(score) <- names(results[[1]][["score"]])

  result <- list(
    projection = projection,
    prediction = prediction,
    confusion = confusion,
    metric = metric,
    probs = probs,
    model = model,
    score = score
  )

  class(result) <- "predict"
  return(result)
}

### Auxiliary functions
# Stack projected blocks as a matrix and set colnames as "block_ncomp"
reformat_projection <- function(projection) {
  names <- names(projection)
  ncomp <- vapply(projection, NCOL, FUN.VALUE = 1L)
  comp_nums <- unlist(lapply(ncomp, seq))
  projection <- as.data.frame(Reduce("cbind", projection))
  colnames(projection) <- paste(
    unlist(Map(function(name, num) rep(name, num), names, ncomp)),
    comp_nums,
    sep = "_"
  )
  return(projection)
}

# Train a model from caret on (X_train, y_train) and make a prediction on
# X_test and evaluate the prediction quality by comparing to y_test.
core_prediction <- function(prediction_model, X_train, X_test,
                            y_train, y_test, metric, classification = FALSE,
                            no_y_test = FALSE, ...) {
  if (classification) {
    y_train <- as.factor(as.matrix(y_train))
    y_test <- factor(as.matrix(y_test), levels = levels(y_train))
  }
  data <- as.data.frame(cbind(X_train, obs = unname(y_train)))
  model <- train(obs ~ .,
    data      = data,
    method    = prediction_model,
    trControl = trainControl(method = "none"),
    na.action = "na.exclude",
    ...
  )

  prediction_train <- data.frame(
    obs = unname(y_train),
    pred = predict(model, X_train)
  )
  prediction_test <- data.frame(
    obs = unname(y_test),
    pred = predict(model, X_test)
  )
  idx_train <- !(is.na(prediction_train$obs) | is.na(prediction_train$pred))
  idx_test <- !(is.na(prediction_test$obs) | is.na(prediction_test$pred))

  probs_train <- probs_test <- NULL

  if (classification) {
    confusion_train <- confusionMatrix(prediction_train$pred,
      reference = prediction_train$obs
    )
    confusion_test <- confusionMatrix(prediction_test$pred,
      reference = prediction_test$obs
    )
    if (is.function(prediction_model$prob)) {
      probs_train <- data.frame(predict(model, X_train, type = "prob"))
      probs_test <- data.frame(predict(model, X_test, type = "prob"))
    }
    metric_train <- multiClassSummary(
      data = prediction_train[idx_train, ],
      lev = levels(prediction_train$obs)
    )
    metric_test <- multiClassSummary(
      data = prediction_test[idx_test, ],
      lev = levels(prediction_test$obs)
    )
  } else {
    confusion_train <- confusion_test <- NA
    metric_train <- postResample(
      pred = prediction_train$pred[idx_train],
      obs = prediction_train$obs[idx_train]
    )
    metric_test <- postResample(
      pred = prediction_test$pred[idx_test],
      obs = prediction_test$obs[idx_test]
    )
  }
  score <- metric_test[grep(metric, names(metric_test), fixed = TRUE)[1]]

  if (no_y_test) {
    score <- confusion_test <- NA
    metric_test <- NULL
    prediction_test$obs <- NULL
  }

  return(list(
    score = score,
    model = model,
    probs = list(train = probs_train, test = probs_test),
    metric = list(train = metric_train, test = metric_test),
    confusion = list(train = confusion_train, test = confusion_test),
    prediction = list(train = prediction_train, test = prediction_test)
  ))
}
