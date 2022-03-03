#' Predict RGCCA
#'
#' Predict a new block from a RGCCA
#'
#' @inheritParams rgcca_transform
#' @param prediction_model A character giving the function used to compare the trained and the tested models
#' @param block_to_predict A character or integer giving the block to predict (must be the same name among train and test set)
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
#' block_to_predict = "industry"
#' y_train = kmeans(A[[block_to_predict]], 3)$cluster
#' y_test = kmeans(X[[block_to_predict]], 3)$cluster
#' res  = rgcca_predict(object, X, block_to_predict = "industry")
#' library(MASS)
#' @importFrom MASS lda
#' @export
rgcca_predict = function(
  rgcca_res,
  X,
  block_to_predict,
  X_scaled = TRUE,
  task = "regression",
  prediction_model = "lm"
  #   regress_on="block" #TODO
) {
  #TODO: Enable to not have the block to predict among the test blocks
  #TODO: Put back asDisjonctive to enable to project block_to_predict without problems
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
  match.arg(task, c("regression", "classification"))
  if (task == "classification" && prediction_model != "lda")
    stop_rgcca("Only \"lda\" model is available for classification task.")
  if (task == "regression" && (prediction_model != "lm" && prediction_model != "cor"))
    stop_rgcca("Only \"cor\" and \"lm\" are available for regression tasks.")
  if (is.null(names(X))) stop_rgcca("Please provide names for the blocks X.")

  ### Check that block_to_predict is among both training and test blocks
  if (is.character(block_to_predict)) {
    train_idx = match(block_to_predict, names(rgcca_res$call$blocks))
    test_idx  = match(block_to_predict, names(X))
  }
  else {
    train_idx = match(block_to_predict, 1:length(rgcca_res$call$blocks))
    test_idx  = match(block_to_predict, 1:length(X))
  }
  if (is.na(train_idx) || is.na(test_idx))
    stop_rgcca(paste0("The block to predict is not among both train and ",
                      "test blocks. Please provide an appropriate one."))
  if (names(X)[[test_idx]] != names(rgcca_res$call$blocks)[[train_idx]])
    stop_rgcca(paste0("Block to predict was provided as an integer but ",
                      "associated block names do not match between train and ",
                      "test blocks. Reorder your blocks or use a name."))

  ### Get train and test target
  y_train = rgcca_res$call$raw[[train_idx]]
  y_test  = X[[test_idx]]

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
      X[[test_idx]] = asDisjonctive(X[[test_idx]], levs = unique(rgcca_res$call$raw[[response]]))
    }
  }

  ### Get projected train and test data
  pred    = rgcca_transform(rgcca_res, X, X_scaled)
  X_train = rgcca_res$Y[names(pred)]
  X_train = reformat_projection(X_train[-test_idx])
  X_test  = reformat_projection(pred[-test_idx])

  # Keep same lines in X_train and y_train
  y_train = subset_rows(y_train, rownames(X_train))

  ### Train prediction model and predict results on X_test
  # TODO: Work on this part
  prediction = NULL
  res        = NULL
  if (task == "regression")
  {

    if (prediction_model == "lm")
    {

      #  if(regress_on=="block")
      #  {
      ychapo <- sapply(
        colnames(y_train),
        function(x) {
          predict(
            lm(
              as.formula(paste(x, " ~ ",paste(colnames(X_train),collapse="+"))),
              data = cbind(X_train, y_train),
              na.action = "na.exclude"
            ),
            cbind(X_test, y_test))
        })
      if (any(is.na(ychapo)))
        warning("NA in predictions.")

      # f defined but never used
      f <- quote(
        if (is.null(dim(y)))
          y[x]
        else
          y[, x]
      )

      prediction=ychapo
      res=data.matrix(y_test - ychapo)

      if (is.null(dim(res))||dim(res)[1]==1) {
        rmse=sqrt(mean(res^2,na.rm=T))
        score <- rmse
      } else{
        rmse<- apply(res,2,function(x){return(sqrt(mean(x^2,na.rm=T)))})
        score <- mean(rmse)
        #score <- mean(apply(res, 2, mean))
      }

      #   }
      #   if(regress_on=="comp")
      #   {
      #TODO
      #   }
    }
    if(prediction_model=="cor")
    {
      # X_test.cor <- get_comp_all(rgcca_res, X=newA, type = "test", pred = pred)
      X_test.cor <- get_comp_all(rgcca_res, X = X, type = "test", pred = pred)

      if (is.null(X[[1]])) {
        # TODO ??? check case for vector
        X_test
      } else{
        connection <- rgcca_res$call$connection[names(pred), names(pred)]
        # cor <- get_cor_all(rgcca_res, newA, X_test)
        # cor <- get_cor_all(rgcca_res, newA3, X_test.cor)
        cor <- get_cor_all(rgcca_res, X, X_test.cor)

        for (i in seq(length(cor))) {
          cor[[i]] <- mean(
            abs(cor[[i]] * connection)[upper.tri(connection)],
            na.rm = TRUE
          )
        }

        score <- mean(unlist(cor), na.rm = TRUE)
      }
    }
    class.fit <- NULL
  }
  if (task == "classification")
  {
    # ngroups defined but never used
    ngroups   <- nlevels(as.factor(y_train))
    class.fit <- switch(prediction_model,
                        "lda"      = {
                          data_for_lda=cbind(X_train,y_train)
                          colnames(data_for_lda)[ncol(data_for_lda)]="quali"
                          reslda     <- lda(quali~., data=data_for_lda, na.action = "na.exclude")
                          class.fit  <- predict(reslda, X_test)$class
                        }#,
                        #             # "logistic" = {
                        #             #     if (ngroups > 2) {
                        #             #         reslog      <- nnet::multinom(y ~ ., data = cbind(X_train, y = y_train), trace = FALSE, na.action = "na.exclude")
                        #             #         class.fit   <- predict(reslog, newdata = cbind(X_test, y = y_test))
                        #             #     } else if (ngroups == 2) {
                        #             #         levs=levels(factor(y_train))
                        #             #         y_train=factor(y_train,levels=levs)
                        #             #         y_test=factor(y_test,levels=levs)
                        #             #         data_for_lda=cbind(X_train,y_train)
                        #             #         colnames(data_for_lda)[ncol(data_for_lda)]="quali"
                        #             #         reslog      <- glm(y ~ ., data = cbind(X_train, y = y_train), family = binomial,na.action="na.exclude")
                        #             #         class.fit   <- predict(reslog, type = "response", newdata = cbind(X_test, y = y_test))
                        #             #         class.fit.class <- class.fit > 0.5 # TODO: cutoff parameter
                        #             #         class.fit       <- factor(class.fit.class)
                        #             #     }
                        #             # }
    )


    if(length(class.fit)==1)
    {
      res=class.fit==as.vector(y_test)
      score=1-res
    }
    else
    {
      res=class.fit ==y_test
      score <- 1-(sum(res) / length(y_test))
    }


  }

  result=list(
    pred = pred,
    #    pred_A=pred_A,
    prediction=prediction,
    class.fit = class.fit,
    score = score,
    res = res,
    rgcca_res=rgcca_res
  )

  class(result)="predict"
  return(result)
}
