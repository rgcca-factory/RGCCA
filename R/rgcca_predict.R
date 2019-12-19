#' Predict RGCCA
#' 
#' Predict a new block from a RGCCA
#' 
#' @inheritParams plot_ind
#' @param newA A list of either a dataframe/matrix or a vector giving the blocks to be predicted
#' @param fit A character giving the function used to compare the trained and the tested models
#' @param bloc_to_pred A character giving the block to predicted (must be the same name among train and test set)
#TODO: either an integer for block_to_pred
#' @param type A character corresponding to the type of prediction among : regression or classification
#' @param y.train A dataframe or a matrix giving the block used as a response in the training
#' @param y.test A dataframe or a matrix giving the block to be predicted
#' @param scale_size_bloc A boolean giving the possibility to scale the blocks by the square root of their column number
#' @param new_scaled A boolean scaling the blocks to predict
#' @param bigA to permeform data reduction for cross-validation, the dataset where A and newA were extracted
#' @examples 
#' library(RGCCA)
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
#' A = lapply(blocks, function(x) x[1:32,])
#' A = lapply(A, function(x) scale2 (x, bias = TRUE) / sqrt(NCOL(x)) )
#' object = rgcca.analyze(A, connection = C, tau = c(0.7,0.8,0.7), 
#'     ncomp = c(3,2,4), superblock = FALSE)
#' newA = lapply(blocks, function(x) x[-c(1:32),])
#' newA = lapply( newA, function(x) x[, sample(1:NCOL(x))] )
#' newA = sample(newA, length(newA))
#' bloc_to_pred = "industry"
#' y.train = kmeans(A[[bloc_to_pred]], 3)$cluster
#' y.test = kmeans(newA[[bloc_to_pred]], 3)$cluster
#' ( res  = rgcca_predict(object, newA, bloc_to_pred = "industry", bigA = blocks) )
#' ( res  = rgcca_predict(object, newA, "regression", "cor", "industry") )
#' ( res  = rgcca_predict(object, newA) )
#' ( res  = rgcca_predict(object, newA = newA, type = "regression", fit = "lm",
#' y.train = A[[bloc_to_pred]], y.test = newA[[bloc_to_pred]] ) )
#' library(MASS)
#' ( res  = rgcca_predict(object, newA = newA, type = "classification",
#' fit = "lda", y.train = y.train, y.test = y.test ) )
#' @importFrom MASS lda
#' @importFrom nnet multinom
#' @export
rgcca_predict = function(
    rgcca,
    newA,
    type = "regression",
    fit = "lm",
    bloc_to_pred = NULL,
    y.train = NULL,
    y.test = NULL,
    bigA = NULL,
    new_scaled = FALSE,
    scale_size_bloc = TRUE) {
    match.arg(type, c("regression", "classification"))
    match.arg(fit, c("lm", "cor", "lda", "logistic"))
    astar <- rgcca$astar
    p <- sapply(rgcca$blocks, ncol)
    B <- length(rgcca$blocks)

    if (type == "classification" && (fit == "cor" || fit == "lm"))
        stop("Please, classification prediction only works with LDA and LOGISTIC")

    if (type == "regression" &&
            (fit == "lda" || fit == "logistic"))
        stop("Please, regression prediction only works with LM and COR")
    
    if (missing(y.train) || missing(y.test)) {
        if (!missing(bloc_to_pred) &&
                !bloc_to_pred %in% names(rgcca$blocks))
            stop("Please, block to predict do not exist")
    }

    # Compute test parameters
    if (is.null(dim(newA[[1]])))
        # case of variable to be predicted
        newp  <- sapply(newA, length)
    else
        # case of blocks to be predicted
        newp <- sapply(newA, ncol)
    newB  <- length(newA)

    # Check similarity between TRAIN and TEST set
    if (is.null(names(rgcca$blocks)) ||  is.null(names(newA)))
        stop("Please, blocs do not have names")

    if (B != newB)
        stop("Please, number of blocs is not the same")

    MATCH <- match(names(newA), names(rgcca$blocks))

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

    MATCH_col <-
        mapply(function(x, y)
            match(get_dim(x)(x), get_dim(y)(y)), newA, rgcca$blocks[MATCH])

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
                x[, y]

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
            
            if (scale_size_bloc)
                res / sqrt(length(data))
            else
                res

        }

        newA <- mapply(
            scl_fun,
            newA,
            reorderList(rgcca$blocks, t_attr = "scaled:center"),
            reorderList(rgcca$blocks, t_attr = "scaled:scale"),
            SIMPLIFY = FALSE
        )

    }


    # Dimension Reduction
    for (i in seq(length(rgcca$blocks)))
        colnames(rgcca$astar[[i]]) <- colnames(rgcca$Y[[i]])
    astar <- reorderList(rgcca$astar, g = TRUE)

    if (is.null(dim(newA[[1]])))
        pred <- lapply(seq(length(newA)), function(x)
            t(as.matrix(newA[[x]])) %*% astar[[x]])
    else
        pred <- lapply(seq(length(newA)), function(x)
            as.matrix(newA[[x]]) %*% astar[[x]])

    if (missing(bloc_to_pred))
        return(list(pred = pred))

    bloc_y <- match(bloc_to_pred, names(rgcca$blocks))

    if (missing(bloc_to_pred) && is.null(colnames(y.train)))
        newbloc_y <- .Machine$integer.max
    else
        newbloc_y <- match(bloc_to_pred, names(newA))


    # Y definition

    if (missing(y.train))
        y.train <- rgcca$blocks[[bloc_y]][, MATCH_col[[newbloc_y]]]

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

    # Prediction
    getComp <- function(type = c("train", "test")) {
        comps <- rgcca$Y[MATCH]
        names <- unlist(lapply(comps[-newbloc_y], colnames))

        if (type ==  "train")
            y <- lapply(seq(length(rgcca$blocks)), function(x)
                comps[[x]][row.names(rgcca$blocks[[x]]), ])
        else
            y <- pred

        # matrix of Y for all selected blocks and all components
        res  <- as.data.frame(Reduce("cbind", y[-newbloc_y]))
        # vector of character giving the name of the block and the number of the component
        col_names <- paste(unlist(mapply(
                function(name, times)  rep(name, times),
                names(newA)[-newbloc_y],
                rgcca$ncomp
            )), names,  sep = "_")
        colnames(res) <- col_names

        return(res)
    }


    rgcca$ncomp <- rgcca$ncomp[MATCH][-newbloc_y]
    comp.train <- getComp("train")
    comp.test <- getComp("test")

    # Scores
    res <- NULL

    if (type == "regression") {
        score <- switch(fit,
            "lm"  = {
                reslm  <- lm(y.train ~ ., data = comp.train, na.action = "na.exclude")
                ychapo <- predict(reslm, comp.test)
                if (any(is.na(ychapo)))
                    warning("NA in predictions.")
                #mean(y.test - ychapo**2)
                # mean(diag(abs(cor(y.test, ychapo))))
                # apply(y.test - ychapo, 2, function(x) mean(abs(x)))
                if (is.null(bigA))
                    bigA <- rgcca$blocks
                m <- function(x)
                    apply(bigA[[bloc_to_pred]], 2, x)

                if (is.null(dim(newA[[1]]))) {
                    f <- quote(y[x])
                } else
                    f <- quote(y[, x])

                r <- function(y)
                        sapply(seq(ncol(bigA[[bloc_to_pred]])), function(x)
                            (eval(f) - m(min)[x]) / (m(max)[x] - m(min)[x]))

                if (is.null(dim(newA[[1]]))) {
                    res <- abs(r(y.test) - r(ychapo))
                    score <- mean(res)
                } else{
                    res <- apply(r(y.test) - r(ychapo), 2, function(x)
                        abs(x))
                    score <- mean(apply(res, 2, mean))
                }


            },
            "cor" = {
                if (is.null(newA[[1]])) {
                    # TODO ??? check case for vector
                    comp.test
                } else{
                    rgcca$C <- rgcca$C[MATCH, MATCH]
                    comp <- list()
                    
                    for (i in seq(max(rgcca$ncomp))) {
                        comp[[i]] <-  matrix(
                            NA,
                            NROW(comp.test),
                            length(newA),
                            dimnames = list(rownames(comp.test), names(newA))
                        )
                        
                        for (n in names(newA)) {
                            pos <- grep(paste0(n, "_comp", i),
                                names(comp.test))
                            if (length(pos) > 0)
                                comp[[i]][, n] <- comp.test[, pos]
                        }
                        comp[[i]] <- sum(abs(
                                cor(comp[[i]], use = "pairwise.complete.obs") * rgcca$C
                            )[upper.tri(rgcca$C)], na.rm = TRUE)
                        if (comp[[i]] == 0)
                            comp[[i]] <-  NA
                        # (cor(comp[[i]], use = "pairwise.complete.obs")*rgcca$C)[upper.tri(rgcca$C)]**2
                    }
                    score <- mean(unlist(comp), na.rm = TRUE)
                }
            })
        class.fit <- NULL
    } else if (type == "classification") {
        ngroups   <- nlevels(as.factor(y.train))
        class.fit <- switch(fit,
            "lda"      = {
                print(comp.train)
                reslda     <- lda(x = comp.train, grouping = y.train)
                class.fit  <- predict(reslda, comp.test)$class
            },
            "logistic" = {
                if (ngroups > 2) {
                    reslog      <- nnet::multinom(y.train ~ ., data = comp.train, trace = FALSE)
                    class.fit   <- predict(reslog, newdata = comp.test)
                } else if (ngroups == 2) {
                    reslog      <- glm(y.train ~ ., data = comp.train, family = binomial)
                    class.fit   <- predict(reslog, type = "response", newdata = comp.test)
                    if (type == "classification") {
                        class.fit.class <- class.fit > 0.5 # TODO: cutoff parameter
                        class.fit       <- factor(as.numeric(class.fit.class))
                    }
                }
            })
        score <- sum(class.fit == y.test) / length(y.test)
    }


    list(
        pred = pred,
        class.fit = class.fit,
        score = score,
        res = res
    )
}
