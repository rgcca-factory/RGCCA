#' Predict RGCCA
#' 
#' Predicts a new block from a RGCCA
#' 
#' @inheritParams plot_ind
#' @param newA A list of either a dataframe/matrix or a vector giving the blocks to be predicted
#' @param fit A character giving the function used to compare the trained and the tested models
#' @param bloc_to_pred A character giving the block to predicted (must be the same name among train and test set)
#TODO: either an integer for block_to_pred
#' @param model A character corresponding to the model of prediction among : regression or classification
#' @param new_scaled A boolean scaling the blocks to predict
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
#' newA = lapply(object1$call$blocks, function(x) x[-c(1:32),])
#' newA = lapply( newA, function(x) x[, sample(1:NCOL(x))] )
#' newA = sample(newA, length(newA))
#' bloc_to_pred = "industry"
#' to_pred_train = kmeans(A[[bloc_to_pred]], 3)$cluster
#' to_pred_test = kmeans(newA[[bloc_to_pred]], 3)$cluster
#' res  = rgcca_predict(object, newA, bloc_to_pred = "industry")
#' # ( res  = rgcca_predict(object, newA, "regression", "cor", "industry") )
#' res  = rgcca_predict(object, newA)
#' library(MASS)
#' @importFrom MASS lda
# @importFrom nnet multinom
#' @export
rgcca_predict = function(
    rgcca_res,
    newA,
    model = "regression",
    fit = "lm",
    bloc_to_pred = NULL,
    new_scaled = TRUE
 #   regress_on="block" #TODO
  ) {
    # Functions to be used
    get_dim <- function(x) {
        if (!is.null(dim(x)))
            return(colnames)
        else
            return(names)
    }
    # Order a a list of matrix or dataframe according to : the index of each element in the list (MATCH);
    # the index of each column in each element (MATCH_COL)
    # l : a list of dataframe or matrix
    # g : a bolean giving the condition to transpose each element of the matrix
    # t_attr : a character giving the name of an attribute of a dataframe or a matrix to order
    reorderList <- function(l, g = FALSE , t_attr = NULL,MATCH,MATCH_col) {
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
    scl_fun <- function(data, center, scale) {
        # Use the scaling parameter of the training set on the tested set
        if(length(center)!=0)
        {
            if (is.null(dim(data)))
                # Case of data is a vector
                data <- t(as.matrix(data))
            res <- scale(data, center, scale) 
        }
        else
        {
            res<-data
        }
        return(res)
    }

    # Checking the input parameters
    if (model == "classification" && (fit == "cor" || fit == "lm"))
        stop_rgcca("Please, classification prediction only works with LDA and LOGISTIC")
    if (model == "regression" &&
        (fit == "lda" || fit == "logistic"))
        stop_rgcca("Please, regression prediction only works with LM and COR")
    stopifnot(is(rgcca_res, "rgcca"))
    match.arg(model, c("regression", "classification"))
    match.arg(fit, c("lm", "cor", "lda", "logistic"))
    for (i in c("new_scaled"))
        check_boolean(i, get(i))
    if (is.null(names(rgcca_res$call$blocks)) ||  is.null(names(newA)))
        stop_rgcca("Please, blocs do not have names")
    
     # Initializations
    prediction=NULL
    if(rgcca_res$call$method=="complete"){rgcca_res$call$blocks=intersection_list(rgcca_res$call$blocks)}
    astar <- rgcca_res$astar
    
    # Dealing with vectors instead of matrices 
    blocks_rgcca_res=lapply(1:length(rgcca_res$call$blocks),
                       function(i) 
                        {
                            if(dim(rgcca_res$call$blocks[[i]])[2]==0 )
                            { 
                                rgcca_res$call$blocks[[i]]=matrix(rgcca_res$call$blocks[[i]],ncol=1);colnames( rgcca_res$call$blocks[[i]])=names(rgcca_res$call$blocks)[i]
                            }
                            if(is.null(colnames(rgcca_res$call$blocks[[i]])))
                            {
                                colnames(rgcca_res$call$blocks[[i]])=names(rgcca_res$call$blocks)[i]
                            }
                            return(rgcca_res$call$blocks[[i]])
                        })
    names(blocks_rgcca_res)=names(rgcca_res$call$blocks)
    newA1=lapply(1:length(newA),
                            function(i) 
                            {
                                if(is.null(dim(newA[[i]])) )
                                { 
                                   newA[[i]]=matrix(newA[[i]],ncol=1);colnames( newA[[i]])=names(newA)[i]
                                }
                                return(newA[[i]])
                            });names(newA1)=names(newA)
    
    # Getting newA2 with disjonctive table (for colnames comparable in response)
    newA2=newA1
    if(!is.null(bloc_to_pred))
    {
        if(mode(newA2[[bloc_to_pred]])=="character")
        {
            if(length(unique(rgcca_res$call$raw[[bloc_to_pred]]))==1){stop("Only one level in the variable to predict")}
            newA2[[bloc_to_pred]]=asDisjonctive(newA2[[bloc_to_pred]],levs=unique(rgcca_res$call$raw[[bloc_to_pred]]))
        }
    }
  
    # Matching the columns in the predict function
     if(model=="regression")
     {  
         MATCH <- match(names(newA1), names(rgcca_res$call$blocks))
         MATCH_col <-  mapply(
             function(x, y) match(get_dim(x)(x), get_dim(y)(y)), 
             newA2, 
             blocks_rgcca_res[MATCH]) 
         MATCH_col2=MATCH_col
     }
     if(model=="classification")
     {
         MATCH <- match(names(newA1), names(rgcca_res$call$raw))
         MATCH_col <-  mapply(
             function(x, y) match(get_dim(x)(x), get_dim(y)(y)), 
             newA1, 
             rgcca_res$call$raw[MATCH])
         MATCH_col2 <-  mapply(
             function(x, y) match(get_dim(x)(x), get_dim(y)(y)), 
             newA2, 
             blocks_rgcca_res[MATCH]) 
     }

    # Checking the matchings
    p <- sapply( blocks_rgcca_res, ncol)
    B <- length(rgcca_res$call$blocks)
    # Matching newA and blocks
    if (sum(is.na(MATCH)) != 0){stop_rgcca("Please, blocs in new data did not exist in old data")}
    
      if (sum(unique(is.na(MATCH_col))) != 0)
        stop_rgcca("Please, some columns names are not the same between the two blocks")
    newp <- sapply(newA2, NCOL)
    newB  <- length(newA2)
    if (sum(newp!=p[MATCH])!=0)
        stop_rgcca("Please, number of columns is not the same")
    if (B != newB)
        stop_rgcca("Please, number of blocks is not the same")
    
    # Scaling de newA (si besoin ie new_scaled=FALSE-cas usuel, et si la rgcca utilise des blocs scales)
    if (!new_scaled  ) { 
        center_vector=reorderList(rgcca_res$call$blocks, t_attr = "scaled:center",MATCH=MATCH,MATCH_col=MATCH_col2)
        scaling_vector=reorderList(rgcca_res$call$blocks, t_attr = "scaled:scale",MATCH=MATCH,MATCH_col=MATCH_col2)
        #center_vector=lapply(rgcca_res$call$blocks,function(x)return(attr(x,"scaled:center")))
        #scaling_vector=lapply(rgcca_res$call$blocks,function(x)return(attr(x,"scaled:scale")))
        # No scaling if  scaling=FALSE, we divide by a vector of ones
        new_scaling_vector=lapply(names(scaling_vector),function(i){
             if(is.null(scaling_vector[[i]]))
            {
                vect_ones=rep(1,length(center_vector[[i]]))
                names(vect_ones)=names(center_vector[[i]])
                return(vect_ones)
            }
            else
            {
                return(scaling_vector[[i]])
            }
            })

        newA3 <- mapply(
            scl_fun,
            newA2,
            center_vector,
            new_scaling_vector,
            SIMPLIFY = FALSE
        )
    } else{ newA3=newA2}
  
# Dimension Reduction
    for (i in seq(length(rgcca_res$call$blocks)))
        colnames(rgcca_res$astar[[i]]) <- colnames(rgcca_res$Y[[i]])
    astar <- reorderList(rgcca_res$astar, g = TRUE,MATCH=MATCH,MATCH_col=MATCH_col2)
    pred <- lapply(seq(length(newA)), function(x)
    {
        M=pm(as.matrix(newA3[[x]]), astar[[x]])
        rownames(M)=rownames(newA[[x]])
        colnames(M)=colnames(astar[[x]])
        return(M)
    }
       )
    names(pred)=names(newA)

    if (missing(bloc_to_pred))
        return(list(pred_y = pred))

    bloc_y <- match(bloc_to_pred, names(rgcca_res$call$blocks))
    newbloc_y <- match(bloc_to_pred, names(newA3))

    # to_pred definition
    if(model=="classification")
    {
        to_pred_train <- rgcca_res$call$raw[[bloc_y]][, MATCH_col[[newbloc_y]], drop = FALSE]
        to_pred_test <- newA[[newbloc_y]]
        
       #  to_pred_train <- rgcca_res$call$raw[[bloc_y]]
        # to_pred_test <- newA[[newbloc_y]]
    }
    if(model=="regression")
    {
       to_pred_train <- blocks_rgcca_res[[bloc_y]][, MATCH_col[[newbloc_y]], drop = FALSE]
        to_pred_test <- newA3[[newbloc_y]] 
       #  to_pred_train <- blocks_rgcca_res[[bloc_y]]
       #  to_pred_test <- newA3[[newbloc_y]] 

        
        if (!is.null(dim(newA[[1]]))) {
            if (any(colnames(to_pred_train) != colnames(to_pred_test)))
                stop_rgcca("Please, train and test sets do not have the same name")
        }
    }
    
    rgcca_res$Y <- rgcca_res$Y[MATCH]
    rgcca_res$call$ncomp <- rgcca_res$call$ncomp[MATCH][-newbloc_y]
    comp.train <- get_comp_all(rgcca_res, newA3, newbloc_y = newbloc_y)
    comp.test <- get_comp_all(rgcca_res, newA3, type = "test", newbloc_y = newbloc_y, pred)

    # Scores
    res <- NULL
    if (model == "regression") 
    {
         
        if(fit=="lm")
        { 
             
              #  if(regress_on=="block")
              #  {
                    ychapo <- sapply(
                        colnames(to_pred_train),
                        function(x) {
                            predict(
                                lm(
                                    as.formula(paste(x, " ~ ",paste(colnames(comp.train),collapse="+"))),
                                    data = cbind(comp.train, to_pred_train), 
                                    na.action = "na.exclude"
                                ), 
                                cbind(comp.test, to_pred_test))
                        })
                    if (any(is.na(ychapo)))
                        warning("NA in predictions.")
                    
                   f <- quote(
                        if (is.null(dim(y)))
                            y[x]
                        else
                            y[, x]
                    )
               
                    prediction=ychapo
                    res=to_pred_test-ychapo
       
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
        if(fit=="cor")
        {
            comp.test.cor <- get_comp_all(rgcca_res, newA=newA, type = "test", pred = pred)
            
                if (is.null(newA3[[1]])) {
                    # TODO ??? check case for vector
                    comp.test
                } else{
                    rgcca_res$call$connection <- rgcca_res$call$connection[MATCH, MATCH]
                   # cor <- get_cor_all(rgcca_res, newA, comp.test)
                    cor <- get_cor_all(rgcca_res, newA3, comp.test.cor)
                    
                    for (i in seq(length(cor))) {
                        cor[[i]] <- mean(
                            abs(
                                cor[[i]] * rgcca_res$call$connection
                            )[upper.tri(rgcca_res$call$connection)], 
                            na.rm = TRUE)     
                    }
                    
                    score <- mean(unlist(cor), na.rm = TRUE)
                }
            }
        class.fit <- NULL
    }
    if (model == "classification")
    {   
        ngroups   <- nlevels(as.factor(to_pred_train))
        class.fit <- switch(fit,
            "lda"      = {
                data_for_lda=cbind(comp.train,to_pred_train)
                colnames(data_for_lda)[ncol(data_for_lda)]="quali"
                reslda     <- lda(quali~., data=data_for_lda, na.action = "na.exclude")
                class.fit  <- predict(reslda, comp.test)$class
            },
            "logistic" = {
                if (ngroups > 2) {
                    reslog      <- nnet::multinom(y ~ ., data = cbind(comp.train, y = to_pred_train), trace = FALSE, na.action = "na.exclude")
                    class.fit   <- predict(reslog, newdata = cbind(comp.test, y = to_pred_test))
                } else if (ngroups == 2) {
                    reslog      <- glm(y ~ ., data = cbind(comp.train, y = to_pred_train), family = binomial,na.action="na.exclude")
                    class.fit   <- predict(reslog, type = "response", newdata = cbind(comp.test, y = to_pred_test))
                    class.fit.class <- class.fit > 0.5 # TODO: cutoff parameter
                    class.fit       <- factor(as.numeric(class.fit.class))
                }
            })
        res=class.fit == to_pred_test
        score <- sum(res) / length(to_pred_test)
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
