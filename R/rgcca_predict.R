#' Predict RGCCA
#'
#' Predict a new block from a fitted RGCCA object.
#'
#' @inheritParams rgcca_transform
#' @param prediction_model A character giving the function used to compare the
#' trained and the tested models.
#' @param response A character or integer giving the block to predict
#' (must be the same name among train and test set).
#' @param ... Additional parameters to be passed to the model fitted on top
#' of RGCCA.
#' @examples
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, 1:3],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' C <- connection <- matrix(
#'   c(
#'     0, 0, 1,
#'     0, 0, 1,
#'     1, 1, 0
#'   ),
#'   3, 3
#' )
#' object1 <- rgcca(blocks,
#'   connection = C, tau = c(0.7, 0.8, 0.7),
#'   ncomp = c(3, 2, 4), superblock = FALSE, response = 3
#' )
#' A <- lapply(object1$call$blocks, function(x) x[1:32, ])
#' object <- rgcca(A,
#'   connection = C, tau = c(0.7, 0.8, 0.7),
#'   ncomp = c(3, 2, 4),
#'   scale = FALSE, scale_block = FALSE,
#'   superblock = FALSE, response = 3
#' )
#' X <- lapply(object1$call$blocks, function(x) x[-c(1:32), ])
#' X <- lapply(X, function(x) x[, sample(1:NCOL(x))])
#' X <- sample(X, length(X))
#' response <- "industry"
#' y_train <- kmeans(A[[response]], 3)$cluster
#' y_test <- kmeans(X[[response]], 3)$cluster
#' res <- rgcca_predict(object, X, response = "industry")
#' @importFrom caret train trainControl confusionMatrix
#' @importFrom caret multiClassSummary postResample
#' @importFrom stats lm predict
#' @export
rgcca_predict <- function(rgcca_res,
                          blocks_test,
                          response,
                          prediction_model = "lm",
                          ...) {
  ### Check input parameters
  if (is.null(names(blocks_test))) {
    stop_rgcca("Please provide names for blocks_test.")
  }

  ### Check that response is among both training and test blocks
  if (is.character(response)) {
    train_idx <- match(response, names(rgcca_res$call$blocks))
    test_idx <- match(response, names(blocks_test))
  } else {
    train_idx <- match(response, seq_along(rgcca_res$call$blocks))
    test_idx <- match(response, seq_along(blocks_test))
  }
  if (is.na(train_idx)) {
    stop_rgcca(paste0(
      "The block to predict is not among train blocks. ",
      "Please provide an appropriate one."
    ))
  }
  if (is.na(test_idx)) {
    no_y_test <- TRUE
    test_idx <- length(blocks_test) + 1
    blocks_test[[names(rgcca_res$call$blocks)[train_idx]]] <-
      rgcca_res$call$raw[[train_idx]]
  } else if (
    names(blocks_test)[[test_idx]] != names(rgcca_res$call$blocks)[[train_idx]]
  ) {
    stop_rgcca(
      "Block to predict was provided as an integer but ",
      "associated block names do not match between train and ",
      "test blocks. Reorder your blocks or use a name."
    )
  } else {
    no_y_test <- FALSE
  }

  tmp <- check_prediction_model(
    prediction_model, rgcca_res$call$raw[[response]]
  )
  prediction_model <- tmp$prediction_model
  classification <- tmp$classification

  ### Get train and test target (if present)
  y_train <- rgcca_res$call$raw[[train_idx]]
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
      stop_rgcca(
        "overfitting model. The RGCCA method led to projected blocks ",
        "that are constant within groups. Try to regularize the model",
        " by increasing tau or decreasing sparsity."
      )
    }
  }

  ### Train prediction model and predict results on X_test
  results <- mapply(
    function(x, y) {
      core_prediction(
        prediction_model, X_train, X_test, x, y, classification, ...
      )
    },
    as.data.frame(y_train),
    y_test
  )

  results <- lapply(
    seq(NCOL(y_train)), function(j) {
      core_prediction(
        prediction_model, X_train, X_test,
        y_train[, j], y_test[, j], classification,
        no_y_test, ...
      )
    }
  )
  names(results) <- colnames(y_train)

  prediction <- as.data.frame(lapply(results, function(res) {
    res[["prediction"]]$test[, "pred"]
  }))

  result <- list(
    projection = projection,
    prediction = prediction,
    results = results,
    rgcca_res = rgcca_res,
    score = mean(unlist(lapply(results, "[[", "score")))
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
    unlist(mapply(function(name, num) rep(name, num), names, ncomp)),
    comp_nums,
    sep = "_"
  )
  return(projection)
}

# Train a model from caret on (X_train, y_train) and make a prediction on
# X_test and evaluate the prediction quality by comparing to y_test.
core_prediction <- function(prediction_model, X_train, X_test,
                            y_train, y_test, classification = FALSE,
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

  if (classification) {
    confusion_train <- confusionMatrix(prediction_train$pred,
      reference = prediction_train$obs
    )
    confusion_test <- confusionMatrix(prediction_test$pred,
      reference = prediction_test$obs
    )
    if (is.null(prediction_model$prob)) {
      prediction_train <- data.frame(cbind(
        prediction_train,
        predict(model, X_train, type = "prob")
      ))
      prediction_test <- data.frame(cbind(
        prediction_test,
        predict(model, X_test, type = "prob")
      ))
    }
    metric_train <- multiClassSummary(
      data = prediction_train,
      lev = levels(prediction_train$obs)
    )
    metric_test <- multiClassSummary(
      data = prediction_test,
      lev = levels(prediction_test$obs)
    )
    score <- 1 - metric_test["Accuracy"]
  } else {
    confusion_train <- confusion_test <- NA
    metric_train <- postResample(
      pred = prediction_train$pred,
      obs = prediction_train$obs
    )
    metric_test <- postResample(
      pred = prediction_test$pred,
      obs = prediction_test$obs
    )
    score <- metric_test["RMSE"]
  }

  if (no_y_test) {
    score <- confusion_test <- NA
    metric_test <- NULL
    prediction_test$obs <- NULL
  }

  return(list(
    model = model,
    prediction = list(train = prediction_train, test = prediction_test),
    metric = list(train = metric_train, test = metric_test),
    confusion = list(train = confusion_train, test = confusion_test),
    score = score
  ))
}
