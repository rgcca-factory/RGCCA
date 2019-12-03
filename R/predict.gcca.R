# object : A list of object giving the rgcca output
# newA : A list of either a dataframe/matrix or a vector giving the blocks to be predicted
# fit : A character giving the function used to compare the trained and the tested models
# block_to_pred : A character giving the block to predicted (must be the same name among train and test set)
#TODO: either an integer for block_to_pred
# y.train : A dataframe or a matrix giving the block used as a response in the training
# y.test : A dataframe or a matrix giving the block to be predicted
# scale_size_bloc : A boolean giving the possibility to scale the blocks by the square root of their column number
# bigA : to permeform data reduction for cross-validation, the dataset where A and newA were extracted
# Examples
# library(RGCCA)
# data("Russett")
# blocks = list(
# agriculture = Russett[, 1:3],
# industry = Russett[, 4:5],
# politic = Russett[, 6:11]
# )
# C = connection = matrix(c(0, 0, 1,
# 0, 0, 1,
# 1, 1, 0),
# 3, 3)
# A = lapply(blocks, function(x) x[1:32,])
# A = lapply(A, function(x) scale2 (x, bias = TRUE) / sqrt(NCOL(x)) )
# object = sgcca(A, C = C, c1 = c(0.7,0.8,0.7), ncomp = c(3,2,4), verbose = F)
# newA = lapply(blocks, function(x) x[-c(1:32),])
# newA = lapply( newA, function(x) x[, sample(1:NCOL(x))] )
# newA = sample(newA, length(newA))
# bloc_to_pred = "industry"
# y.train = kmeans(A[[bloc_to_pred]], 3)$cluster
# y.test = kmeans(newA[[bloc_to_pred]], 3)$cluster
# ( res  = predict.gcca(object, A, newA, "regression", "lm", "industry", bigA = blocks) )
# ( res  = predict.gcca(object, A, newA, "regression", "cor", "industry") )
# ( res  = predict.gcca(object, A ) )
# ( res  = predict.gcca(object, A, newA = newA, type = "regression", fit = "lm", y.train = A[[bloc_to_pred]], y.test = newA[[bloc_to_pred]] ) )
# ( res  = predict.gcca(object, A, newA = newA, type = "classification", fit = "lda", y.train = y.train, y.test = y.test ) )
#' @importFrom MASS lda
#' @importFrom nnet multinom


predict.gcca = function(
    object,
    A,
    newA,
    type = c("regression", "classification"),
    fit = c("lm", "cor", "lda", "logistic"),
    bloc_to_pred = NULL,
    y.train = NULL,
    y.test = NULL,
    bigA = NULL,
    new_scaled = FALSE,
    scale_size_bloc = TRUE,
    cutoff = 0.5) {

  type = match.arg(type)
  astar = object$astar
  p = sapply(A, ncol)
  B = length(A)

      # Arguments checking
      if (missing(fit) || missing(type) )
        stop("Please, define type, fit and bloc_to_pred or y.train/y.test")

      if (type == "classification" && (fit == "cor" || fit == "lm"))
        stop("Please, classification prediction only works with LDA and LOGISTIC")

      if (type == "regression" && (fit == "lda" || fit == "logistic"))
        stop("Please, regression prediction only works with LM and COR")

      if ( missing(bloc_to_pred) && ( missing(y.train) || missing(y.test) ) )
        stop("Please, fill either bloc_to_pred or y.train and y.test parameters")

      # Compute test parameters
      if (is.null(dim(newA[[1]])))
        # case of variable to be predicted
        newp  = sapply(newA, length)
      else
        # case of blocks to be predicted
        newp  = sapply(newA, ncol)
      newB  = length(newA)

      # Check similarity between TRAIN and TEST set
      if (is.null(names(A)) ||  is.null(names(newA)))
        stop("Please, blocs do not have names")

      if (length(A) != B && length(newA) != newB)
        stop("Please, number of blocs is not the same")

      MATCH = match(names(newA), names(A))

      if (sum(is.na(MATCH)) != 0)
        stop("Please, blocs in new data did not exist in old data")

      if (!identical(newp, p[MATCH]))
        stop("Please, number of column is not the same")
      
      if (!bloc_to_pred %in% names(A))
         stop("Please, block to predict do not exist")


      # Y definition
      if (missing(bloc_to_pred)) {
          bloc_to_pred = names (which (unlist (lapply(A, function(i) {
              all(unlist(lapply(colnames(y.train),
                  function(j)
                      j %in% colnames(i))))
          }))))
      }

      bloc_y = match(bloc_to_pred, names(A))
      newbloc_y = match(bloc_to_pred, names(newA))
      
      f1 <- f2 <- names
      if (!is.null(dim(newA[[1]])))
          f1 <- colnames
      if (!is.null(dim(A[[1]])))
          f2 <- colnames
      
      MATCH_col = mapply(function(x, y) match(f1(x), f2(y)), newA, A[MATCH])
      
      if (sum(unique(is.na(MATCH_col))) != 0)
        stop("Please, some columns names are not the same between the two blocks")


      # Order a a list of matrix or dataframe according to : the index of each element in the list (MATCH);
      # the index of each column in each element (MATCH_COL)
      # l : a list of dataframe or matrix
      # g : a bolean giving the condition to transpose each element of the matrix
      # t_attr : a character giving the name of an attribute of a dataframe or a matrix to order
      reorderList = function(l, g = FALSE , t_attr = NULL){

        # Deals with a transposed matrix
        if (!is.null(t_attr))
          f_attr = function(x, y) attr(x, t_attr)[y]
        else
          f_attr = function(x, y) x[, y]

        # Deals with attributes from a dataframe
        if (isTRUE(g))
          g = function(x) t(x)
        else
          g = function(x) x

        mapply(function(x, y) g(f_attr(g(x), y)), l[MATCH], MATCH_col)
      }

      # TODO: if new as an only individuals
      # scale before y_train and y_test attribution ? Otherwise, the response is not scaled
      # Scaling
      if (!new_scaled) {

        scl_fun = function(data, center, scale) {
          # Use the scaling parameter of the training set on the tested set

          if (is.null(dim(data)))
            # Case of data is a vector
            data = t(as.matrix(data))

          res = scale(data, center, scale)

          if (scale_size_bloc)
            res / sqrt(length(data))
          else
            res

        }

        newA = mapply(scl_fun,
                      newA,
                      reorderList(A, t_attr = "scaled:center"),
                      reorderList(A, t_attr = "scaled:scale"),
                      SIMPLIFY = FALSE)

      }

      # Y definition

      if (missing(y.train))
        y.train = A[[bloc_y]][, MATCH_col[[newbloc_y]] ]

      # TODO : sampled columns for y.test

      if (missing(y.test)) {
        # if (length(newbloc_y) < 1)
        #   warning("Please, bloc_to_pred and y.test are absent from testing set. Coeffcients only will be given in outpus.") #TODO
        y.test = newA[[newbloc_y]]
      }else{
        # if (is.na(newbloc_y))
        #   warning("Please, bloc_to_pred is absent from testing set")
        MATCH_col_y = match(colnames(y.train), colnames(y.test))
        y.test = y.test[MATCH_col_y]
      }

      if (!is.null(dim(newA[[1]]))) {
        if ( any(colnames(y.train) != colnames(y.test) ) )
          stop("Please, train and test sets do not have the same name")
      }

      # Dimension Reduction
      for (i in 1:length(A))
        colnames(object$astar[[i]]) = colnames(object$Y[[i]])
      astar = reorderList(object$astar, g = TRUE)

        if (is.null(dim(newA[[1]])))
        pred = lapply(1:length(newA), function(x) t(as.matrix(newA[[x]])) %*% astar[[x]])
      else
        pred = lapply(1:length(newA), function(x) as.matrix(newA[[x]]) %*% astar[[x]])


      # Prediction

      getComp = function(type = c("train", "test")){

          comps <- object$Y[MATCH]
          names <- unlist(lapply(comps[-newbloc_y], colnames))
          
        if (type ==  "train")
          y = lapply(1:length(A), function(x) comps[[x]][row.names(A[[x]]),])
        else
          y = pred

        # matrix of Y for all selected blocks and all components
        res  = as.data.frame( Reduce("cbind", y[-newbloc_y]) )
        # vector of character giving the name of the block and the number of the component
        col_names = paste( unlist(mapply(function(name, times) rep(name, times),
                                         names(newA)[-newbloc_y],
                                         object$ncomp)),
                           names,
                           sep = "_")
        colnames(res) = col_names

        return(res)
      }

    object$ncomp = object$ncomp[MATCH][-newbloc_y]
    comp.train = getComp("train")
    comp.test = getComp("test")

      # Scores

      res <- NULL
    
      if (type == "regression") {
        score = switch(
          fit,
          "lm"  = {
            reslm  = lm(y.train ~ ., data = comp.train, na.action = "na.exclude")
            ychapo = predict(reslm, comp.test)
            if (any(is.na(ychapo))) warning("NA in predictions.")
            #mean(y.test - ychapo**2)
            # mean(diag(abs(cor(y.test, ychapo))))
            # apply(y.test - ychapo, 2, function(x) mean(abs(x)))
            if (is.null(bigA))
                bigA <- A
            m <- function(x) apply(bigA[[bloc_to_pred]], 2, x)
            
            if (is.null(dim(newA[[1]]))){
                f <- quote(y[x])
            }else
                f <- quote(y[, x])
            
            r <- function(y) sapply(1:ncol(bigA[[bloc_to_pred]]), function(x) (eval(f) - m(min)[x]) / (m(max)[x] - m(min)[x]))
            
            if (is.null(dim(newA[[1]]))){
                res <- abs(r(y.test) - r(ychapo))
                 score <- mean(res)
            }else{
                res <- apply(r(y.test) - r(ychapo), 2, function(x) abs(x))
                score <- mean(apply(res, 2, mean))
            }
            
            
          },
          "cor" = {
            if (is.null(newA[[1]])) { # TODO ??? check case for vector
              comp.test
            }else{
              object$C = object$C[MATCH, MATCH]
              comp = list()
              
              for (i in 1:max(object$ncomp)) {

                comp[[i]] =  matrix(NA,
                                    NROW(comp.test),
                                    length(newA),
                                    dimnames = list(rownames(comp.test), names(newA)))

                for (n in names(newA)) {
                  pos = grep(paste0(n, "_comp", i), names(comp.test) )
                  if (length(pos) > 0)
                    comp[[i]][, n] = comp.test[, pos]
                }
                comp[[i]] = sum(abs(cor(comp[[i]], use = "pairwise.complete.obs")*object$C)[upper.tri(object$C)], na.rm = TRUE)
                if (comp[[i]] == 0)
                    comp[[i]] =  NA
                # (cor(comp[[i]], use = "pairwise.complete.obs")*object$C)[upper.tri(object$C)]**2
              }
              score <- mean(unlist(comp), na.rm = TRUE)
            }
          }
        )
        class.fit = NULL
      } else if (type == "classification") {
        ngroups   = nlevels(as.factor(y.train))
        class.fit = switch(
          # K = 3
          # y.train = kmeans(y.train, K)$cluster
          # y.test = kmeans(y.test, K)$cluster
          fit,
          "lda"      = {
           # library("MASS")
            reslda     = lda(x = comp.train, grouping = y.train)
            class.fit  = predict(reslda, comp.test)$class
          },
          "logistic" = {
            # K = 2
            # y.train = kmeans(y.train, K)$cluster -1
            # y.test =  kmeans(y.test, K)$cluster -1
            if (ngroups > 2) {
             # library("nnet")
              reslog      = nnet::multinom(y.train ~ ., data = comp.train, trace = FALSE)
              class.fit   = predict(reslog, newdata = comp.test)
            }else if (ngroups == 2) {
              reslog      = glm(y.train ~ ., data = comp.train, family = binomial)
              class.fit   = predict(reslog, type = "response", newdata = comp.test)
              if (type == "classification") {
                class.fit.class = class.fit > cutoff # TODO ???
                class.fit       = factor(as.numeric(class.fit.class))
              }
            }
          }
        )
        score = sum(class.fit == y.test)/length(y.test)
      }
 

  list(pred = pred, class.fit = class.fit, score = score, res = res)
}
