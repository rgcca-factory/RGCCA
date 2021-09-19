#' Predict RGCCA
#'
#' Predict a new block from a RGCCA
#'
#' @inheritParams rgcca_transform
#' @param prediction_model A character giving the function used to compare the trained and the tested models
#' @param response A character or integer giving the block to predict (must be the same name among train and test set)
#' @param task A character corresponding to the prediction task among : regression or classification
#' @examples
#' data("Russett")
#' blocks = list(
#' agriculture = Russett[, 1:3],
#' industry = Russett[, 4:5],
#' politic = Russett[, 6:11]
#' )
#' C = connection = matrix(c(0, 0, 1,
#' 0, 0, 1,
#' 1, 1, 0),
#' 3, 3)
#' object1 = rgcca(blocks, connection = C, tau = c(0.7,0.8,0.7),
#'     ncomp = c(3,2,4), superblock = FALSE, response = 3)
#' A = lapply(object1$call$blocks, function(x) x[1:32,])
#' object = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
#'     ncomp = c(3,2,4), scale = FALSE, scale_block = FALSE, superblock = FALSE, response = 3)
#' X = lapply(object1$call$blocks, function(x) x[-c(1:32),])
#' X = lapply( X, function(x) x[, sample(1:NCOL(x))] )
#' X = sample(X, length(X))
#' response = "industry"
#' y_train = kmeans(A[[response]], 3)$cluster
#' y_test = kmeans(X[[response]], 3)$cluster
#' res  = rgcca_predict(object, X, response = "industry")
#' library(MASS)
#' @importFrom MASS lda
# @importFrom nnet multinom
#' @export
rgcca_predict = function(
  rgcca_res,
  blocks_test,
  response,
  method = "lm",
  ...
  #   regress_on="block" #TODO
) {
  #TODO: Enable to not have the block to predict among the test blocks
  #TODO: Put back asDisjonctive to enable to project response without problems
  ### Auxiliary function
  # Stack projected blocks as a matrix and set colnames as "block_ncomp"
  reformat_projection = function(projection) {
    names      = names(projection)
    ncomp      = sapply(projection, ncol)
    comp_nums  = unlist(sapply(ncomp, seq))
    projection = as.data.frame(Reduce("cbind", projection))
    colnames(projection) <- paste(
      unlist(mapply(function(name, num)  rep(name, num), names, ncomp)),
      comp_nums,
      sep = "_")
    return(projection)
  }

  ### Check input parameters
  all_methods = caret::modelLookup()
  match.arg(method, c("cor", all_methods$model))
  if (is.null(names(blocks_test))) stop_rgcca("Please provide names for the blocks_test.")

  ### Check that response is among both training and test blocks
  if (is.character(response)) {
    train_idx = match(response, names(rgcca_res$call$blocks))
    test_idx  = match(response, names(blocks_test))
  }
  else {
    train_idx = match(response, 1:length(rgcca_res$call$blocks))
    test_idx  = match(response, 1:length(blocks_test))
  }
  if (is.na(train_idx) || is.na(test_idx))
    stop_rgcca(paste0("The block to predict is not among both train and ",
                      "test blocks. Please provide an appropriate one."))
  if (names(blocks_test)[[test_idx]] != names(rgcca_res$call$blocks)[[train_idx]])
    stop_rgcca(paste0("Block to predict was provided as an integer but ",
                      "associated block names do not match between train and ",
                      "test blocks. Reorder your blocks or use a name."))

  ### Get train and test target
  y_train = rgcca_res$call$raw[[train_idx]]
  y_test  = blocks_test[[test_idx]]

  if ((is.null(dim(y_test)) && !all(dim(y_train)[-1] == 1)) ||
      any(dim(y_test)[-1] != dim(y_train)[-1])
  ) stop_rgcca(paste0("Dimensions of block to predict do not match between",
                      " train and test blocks."))
  if (!is.null(dim(y_test))) y_test = y_test[, colnames(y_train), drop = FALSE]

  ### One hot encode block to predict if needed
  response = rgcca_res$call$response
  if (!is.null(response)) {
    if (mode(rgcca_res$call$raw[[response]]) == "character") {
      if (length(unique(rgcca_res$call$raw[[response]])) == 1) {
        stop("Only one level in the variable to predict")}
      blocks_test[[test_idx]] = asDisjonctive(blocks_test[[test_idx]], levs = unique(rgcca_res$call$raw[[response]]))
    }
  }

  ### Get projected train and test data
  blocks_test_scaled = !(rgcca_res$call$scale)
  projection         = rgcca_transform(rgcca_res, blocks_test, blocks_test_scaled)
  X_train            = rgcca_res$Y[names(projection)]
  X_train            = reformat_projection(X_train[-test_idx])
  X_test             = reformat_projection(projection[-test_idx])

  # Keep same lines in X_train and y_train
  y_train = subset_rows(y_train, rownames(X_train))

  ### Train prediction model and predict results on X_test
  core_prediction = function(model_info, X_train, X_test, y_train, y_test){
    data = as.data.frame(cbind(X_train, obs = unname(y_train)))
    model = train(obs   ~ .,
                  data      = data,
                  method    = method,
                  trControl = trainControl(method = "none"), ...)

    prediction_train = data.frame(obs  = unname(y_train),
                                  pred = predict(model, X_train))
    prediction_test  = data.frame(obs  = unname(y_test),
                                  pred = predict(model, X_test))
    if ((model_info$forClass) && (!is.numeric(y_train))){
      prediction_train$obs = as.factor(prediction_train$obs)
      prediction_test$obs  = as.factor(prediction_test$obs)
      confusion_train      = confusionMatrix(prediction_train$pred,
                                             reference = prediction_train$obs)
      confusion_test       = confusionMatrix(prediction_test$pred,
                                             reference = prediction_test$obs)
      if (model_info$probModel){
        prediction_train = data.frame(cbind(prediction_train,
                                            predict(model, X_train, type = "prob")))
        prediction_test  = data.frame(cbind(prediction_test,
                                            predict(model, X_test, type = "prob")))
      }
      metric_train     = multiClassSummary(data = prediction_train,
                                           lev  = levels(prediction_train$obs))
      metric_test      = multiClassSummary(data = prediction_test,
                                           lev  = levels(prediction_test$obs))
    }

    if ((model_info$forReg) && (is.numeric(y_train))){
      confusion_train = confusion_test = NA
      metric_train    = postResample(pred = prediction_train$pred,
                                     obs  = prediction_train$obs)
      metric_test     = postResample(pred = prediction_test$pred,
                                     obs  = prediction_test$obs)
    }
    return(list(model      = model,
                prediction = list(train = prediction_train, test = prediction_test),
                metric     = list(train = metric_train, test = metric_test),
                confusion  = list(train = confusion_train, test = confusion_test))
    )
  }

  if (method != "cor"){
    model_info = all_methods[which(all_methods$model == method), ]
    if ((dim(y_train)[2] >= 2) && (sum(!apply(y_train, 2, is.numeric)) != 0)){
      stop_rgcca(paste0("Multiple columns prediction is not handled."))
    }else{
      results = mapply(function(x, y) core_prediction(model_info, X_train, X_test, x, y),
                       as.data.frame(y_train),
                       as.data.frame(y_test))
    }
  }else if(method == "cor"){
    X_train_all = reformat_projection(rgcca_res$Y[names(projection)])
    X_test_all  = reformat_projection(projection)
    connection  = rgcca_res$call$connection[names(projection), names(projection)]
    connection[connection == 0] = NA
    cor_train   = get_cor_all(rgcca_res = rgcca_res, comps = X_train_all)
    cor_test    = get_cor_all(rgcca_res = rgcca_res, blocks = blocks_test,
                              comps = X_test_all)

    for (i in seq(length(cor))) {
      cor_train[[i]] = mean(abs(cor_train[[i]] * connection)[upper.tri(connection)],
                            na.rm = TRUE)
      cor_test[[i]]  = mean(abs(cor_test[[i]] * connection)[upper.tri(connection)],
                            na.rm = TRUE)
    }
    metric_train = mean(unlist(cor_train), na.rm = TRUE)
    metric_test  = mean(unlist(cor_test), na.rm = TRUE)

    results           = matrix(list(NA), 4, 1)
    results[3, 1]     = list(list(train = metric_train, test = metric_test))
    rownames(results) = c("model", "prediction", "metric", "confusion")
  }

  result=list(projection = projection,
              results    = results,
              rgcca_res  = rgcca_res
              )

  class(result)="predict"
  return(result)
}
