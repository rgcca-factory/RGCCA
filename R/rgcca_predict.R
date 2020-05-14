#' Predict RGCCA
#' 
#' Predict a new block from a RGCCA
#' 
#' @inheritParams plot_ind
#' @param newA A list of either a dataframe/matrix or a vector giving the blocks to be predicted
#' @param fit A character giving the function used to compare the trained and the tested models
#' @param bloc_to_pred A character giving the block to predicted (must be the same name among train and test set)
#TODO: either an integer for block_to_pred
#' @param model A character corresponding to the model of prediction among : regression or classification
#' @param y.train A dataframe or a matrix giving the block used as a response in the training
#' @param y.test A dataframe or a matrix giving the block to be predicted
#' @param scale_size_bloc A boolean giving the possibility to scale the blocks by the square root of their column number
#' @param new_scaled A boolean scaling the blocks to predict
#' @param bigA to permeform data reduction for cross-validation, the dataset where A and newA were extracted
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
#'     ncomp = c(3,2,4), scale = FALSE, sameBlockWeight = FALSE, superblock = FALSE, response = 3)
#' newA = lapply(object1$call$blocks, function(x) x[-c(1:32),])
#' newA = lapply( newA, function(x) x[, sample(1:NCOL(x))] )
#' newA = sample(newA, length(newA))
#' bloc_to_pred = "industry"
#' y.train = kmeans(A[[bloc_to_pred]], 3)$cluster
#' y.test = kmeans(newA[[bloc_to_pred]], 3)$cluster
#' ( res  = rgcca_predict(object, newA, bloc_to_pred = "industry", bigA = blocks) )
#' ( res  = rgcca_predict(object, newA, "regression", "cor", "industry") )
#' ( res  = rgcca_predict(object, newA) )
#' ( res  = rgcca_predict(object, newA = newA, model = "regression", fit = "lm",
#' y.train = A[[bloc_to_pred]], y.test = newA[[bloc_to_pred]] ) )
#' library(MASS)
#' ( res  = rgcca_predict(object, newA = newA, model = "classification",
#' fit = "lda", y.train = y.train, y.test = y.test ) )
#' @importFrom MASS lda
# @importFrom nnet multinom
#' @export
rgcca_predict = function(
    rgcca_res,
    newA,
    model = "regression",
    fit = "lm",
    bloc_to_pred = NULL,
    y.train = NULL,
    y.test = NULL,
    bigA = NULL,
    new_scaled = TRUE,
    scale_size_bloc = TRUE) {

    stopifnot(is(rgcca_res, "rgcca"))
    match.arg(model, c("regression", "classification"))
    match.arg(fit, c("lm", "cor", "lda", "logistic"))
    # if(!is.null(bigA))
    #     check_blocks(bigA)

    for (i in c("new_scaled", "scale_size_bloc"))
            check_boolean(i, get(i))

    astar <- rgcca_res$astar
    p <- sapply(rgcca_res$call$blocks, ncol)
    B <- length(rgcca_res$call$blocks)

    if (model == "classification" && (fit == "cor" || fit == "lm"))
        stop("Please, classification prediction only works with LDA and LOGISTIC")

    if (model == "regression" &&
            (fit == "lda" || fit == "logistic"))
        stop("Please, regression prediction only works with LM and COR")

    if (missing(y.train) || missing(y.test)) {
        if (!missing(bloc_to_pred) &&
                !bloc_to_pred %in% names(rgcca_res$call$blocks))
            stop("Please, block to predict do not exist")
    }

    # Compute test parameters
    if (is.null(dim(newA[[1]])))
        # case of variable to be predicted
        newp  <- sapply(newA, length)
    else
        # case of blocks to be predicted
        newp <- sapply(newA, NCOL)
    newB  <- length(newA)

    # Check similarity between TRAIN and TEST set
    if (is.null(names(rgcca_res$call$blocks)) ||  is.null(names(newA)))
        stop("Please, blocs do not have names")

    if (B != newB)
        stop("Please, number of blocs is not the same")

    MATCH <- match(names(newA), names(rgcca_res$call$blocks))

    if (sum(is.na(MATCH)) != 0)
        stop("Please, blocs in new data did not exist in old data")

    if (!identical(newp, p[MATCH]))
        stop("Please, number of column is not the same")

    get_dim <- function(x) {
        if (!is.null(dim(x)))
            return(colnames)
        else
            return(names)
    }

    MATCH_col <-  mapply(
        function(x, y) match(get_dim(x)(x), get_dim(y)(y)), 
        newA, 
        rgcca_res$call$blocks[MATCH])

    if (sum(unique(is.na(MATCH_col))) != 0)
        stop("Please, some columns names are not the same between the two blocks")

    # Order a a list of matrix or dataframe according to : the index of each element in the list (MATCH);
    # the index of each column in each element (MATCH_COL)
    # l : a list of dataframe or matrix
    # g : a bolean giving the condition to transpose each element of the matrix
    # t_attr : a character giving the name of an attribute of a dataframe or a matrix to order
    reorderList <- function(l, g = FALSE , t_attr = NULL) {
        # Deals with a transposed matrix
        if (!is.null(t_attr))
            f_attr <- function(x, y)
                attr(x, t_attr)[y]
        else
            f_attr <- function(x, y)
                x[, y, drop = FALSE]

        # Deals with attributes from a dataframe
        if (isTRUE(g))
            g <- function(x)
                t(x)
        else
            g <- function(x)
                x

        mapply(function(x, y)
            g(f_attr(g(x), y)), l[MATCH], MATCH_col)
    }

    # TODO: if new as an only individuals
    # scale before y_train and y_test attribution ? Otherwise, the response is not scaled
    # Scaling
    if (!new_scaled) {
        scl_fun <- function(data, center, scale) {
            # Use the scaling parameter of the training set on the tested set
            
            if (is.null(dim(data)))
                # Case of data is a vector
                data <- t(as.matrix(data))
            
            res <- scale(data, center, scale)
            
            # if (scale_size_bloc)
            #     res / sqrt(length(data))
            # else
            #     res

        }

        newA <- mapply(
            scl_fun,
            newA,
            reorderList(rgcca_res$call$blocks, t_attr = "scaled:center"),
            reorderList(rgcca_res$call$blocks, t_attr = "scaled:scale"),
            SIMPLIFY = FALSE
        )

    }


    # Dimension Reduction
    for (i in seq(length(rgcca_res$call$blocks)))
        colnames(rgcca_res$astar[[i]]) <- colnames(rgcca_res$Y[[i]])
    astar <- reorderList(rgcca_res$astar, g = TRUE)

    pred <- lapply(seq(length(newA)), function(x)
    {
        M=as.matrix(newA[[x]]) %*% astar[[x]]
        rownames(M)=rownames(newA[[x]])
        colnames(M)=colnames(astar[[x]])
        return(M)
    }
    )
    names(pred)=names(newA)

    if (missing(bloc_to_pred))
        return(list(pred = pred))

    bloc_y <- match(bloc_to_pred, names(rgcca_res$call$blocks))

    if (missing(bloc_to_pred) && is.null(colnames(y.train)))
        newbloc_y <- .Machine$integer.max
    else
        newbloc_y <- match(bloc_to_pred, names(newA))


    # Y definition

    if (missing(y.train))
        y.train <- rgcca_res$A[[bloc_y]][, MATCH_col[[newbloc_y]], drop = FALSE]

    # TODO : sampled columns for y.test

    if (missing(y.test)) {
        # if (length(newbloc_y) < 1)
        #   warning("Please, bloc_to_pred and y.test are absent from testing set. Coeffcients only will be given in outpus.") #TODO
        y.test <- newA[[newbloc_y]]
    } else{
        # if (is.na(newbloc_y))
        #   warning("Please, bloc_to_pred is absent from testing set")
        # MATCH_col_y <- match(get_dim(y.train)(y.train), get_dim(y.test)(y.test))
        # y.test <- y.test[MATCH_col_y]
    }

    if (!is.null(dim(newA[[1]]))) {
        if (any(colnames(y.train) != colnames(y.test)))
            stop("Please, train and test sets do not have the same name")
    }

    rgcca_res$Y <- rgcca_res$Y[MATCH]
    rgcca_res$call$ncomp <- rgcca_res$call$ncomp[MATCH][-newbloc_y]
    comp.train <- get_comp_all(rgcca_res, newA, newbloc_y = newbloc_y)
    comp.test <- get_comp_all(rgcca_res, newA, type = "test", newbloc_y = newbloc_y, pred)

    # Scores
    res <- NULL

    if (model == "regression") {
        score <- switch(fit,
            "lm"  = {

                ychapo <- sapply(
                    colnames(y.train),
                    function(x) {
                        predict(
                            lm(
                                as.formula(paste(x, " ~ ",paste(colnames(comp.train),collapse="+"))),
                                data = cbind(comp.train, y.train), 
                                na.action = "na.exclude"
                            ), 
                            cbind(comp.test, y.test))
                    })
                
                if (any(is.na(ychapo)))
                    warning("NA in predictions.")
                
                #   if (is.null(bigA))
                #       bigA <- rgcca_res$call$blocks
                
                #   m <- function(x)
                #       apply(bigA[[bloc_to_pred]], 2, x) # TODO
                
                f <- quote(
                    if (is.null(dim(y)))
                        y[x]
                    else
                        y[, x]
                )
                
                #   r <- function(y)
                #       sapply(seq(NCOL(bigA[[bloc_to_pred]])), function(x)
                #           (eval(f) - m(min)[x]) / (m(max)[x] - m(min)[x]))
                
                #   res <- r(y.test) - r(ychapo)
                prediction=ychapo
                res=y.test-ychapo
                if (is.null(dim(res))) {
                    #res <- abs(res)
                    score <- sqrt(mean(res^2))
                } else{
                    #  res <- apply(res, 2, function(x) abs(x))
                    res2 <- apply(res,2,function(x){return(sqrt(mean(x^2)))})
                    score <- sqrt(mean(res2^2))
                    #score <- mean(apply(res, 2, mean))
                }


            },
            "cor" = {
                if (is.null(newA[[1]])) {
                    # TODO ??? check case for vector
                    comp.test
                } else{
                    rgcca_res$call$connection <- rgcca_res$call$connection[MATCH, MATCH]
                    cor <- get_cor_all(rgcca_res, newA, comp.test)

                    for (i in seq(length(cor))) {
                        cor[[i]] <- mean(
                            abs(
                                cor[[i]] * rgcca_res$call$connection
                            )[upper.tri(rgcca_res$call$connection)], 
                            na.rm = TRUE)     
                    }
                    score <- mean(unlist(cor), na.rm = TRUE)
                }
            })
        class.fit <- NULL

    } else if (model == "classification") {
        ngroups   <- nlevels(as.factor(y.train))
        class.fit <- switch(fit,
            "lda"      = {
                reslda     <- lda(x = comp.train, grouping = y.train, na.action = "na.exclude")
                class.fit  <- predict(reslda, comp.test)$class
            },
            "logistic" = {
                if (ngroups > 2) {
                    reslog      <- nnet::multinom(y ~ ., data = cbind(comp.train, y = y.train), trace = FALSE, na.action = "na.exclude")
                    class.fit   <- predict(reslog, newdata = cbind(comp.test, y = y.test))
                } else if (ngroups == 2) {
                    reslog      <- glm(y ~ ., data = cbind(comp.train, y = y.train), family = binomial)
                    class.fit   <- predict(reslog, type = "response", newdata = cbind(comp.test, y = y.test))
                    class.fit.class <- class.fit > 0.5 # TODO: cutoff parameter
                    class.fit       <- factor(as.numeric(class.fit.class))
                }
            })
        score <- sum(class.fit == y.test) / length(y.test)
    }


    result=list(
        pred = pred,
        class.fit = class.fit,
        score = score,
        res = res,
        rgcca_res=rgcca_res
    )
    class(result)="predict"
    return(result)
}
