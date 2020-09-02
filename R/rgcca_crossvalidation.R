#' Cross-validation
#' 
#' Uses cross-validation to validate a predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams rgcca
#' @inheritParams bootstrap
#' @inheritParams plot_ind
#' @param k when k fold is chosen, the k parameter.
#' @param validation Among "loo", "kfold", "test".
#' @param parallelization if TRUE parallelization is run, if FALSE, no parallelisation is run. If NULL (default) parallelization is always used except for Windows in case of length(nperm)<10
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, response = 3,superblock=FALSE)
#' res=rgcca_crossvalidation(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
#' rgcca_crossvalidation(rgcca_out, n_cores = 1)
#' @export
#' @seealso \link{rgcca}, \link{rgcca_predict}, \link{plot.predict}
rgcca_crossvalidation <- function(
    rgcca_res,
    validation = "kfold",
    model = "regression",
    fit = "lm",
    new_scaled = TRUE,
    k = 5,
    scale=NULL,
    scale_block=NULL,
    tol=1e-8,
    scheme=NULL,
    method=NULL,
    type=NULL,
    init=NULL,
    bias=NULL,
    connection=NULL,
    ncomp=NULL,
    tau=NULL,
    sparsity=NULL,
    n_cores = parallel::detectCores() - 1,
    parallelization=TRUE,
    ...) {

    if(is.null(connection)){    connection <- rgcca_res$call$connection }
    if(is.null(scale)){    scale <- rgcca_res$call$scale }
    if(is.null(scale_block)){    scale_block <- rgcca_res$call$scale_block }
    if(is.null(method)){    method <- rgcca_res$call$method }
    if(is.null(scheme)){     scheme <- rgcca_res$call$scheme}
    if(is.null(bias)){     bias <- rgcca_res$call$bias}
    if(is.null(type)){   type <- rgcca_res$call$type}
    if(is.null(init)){     init <- rgcca_res$call$init}
    if(is.null(ncomp)){        ncomp <- rgcca_res$call$ncomp}
    if(is.null(tau)){        tau <- rgcca_res$call$tau}
    if(is.null(sparsity)){        sparsity<- rgcca_res$call$sparsity}
    
    stopifnot(is(rgcca_res, "rgcca"))
    if(is.null(rgcca_res$call$response)){
       stop_rgcca("This function required an analysis in a supervised mode")}
   
    match.arg(validation, c("loo", "test", "kfold"))
    check_integer("k", k, min = 2)
    check_integer("n_cores", n_cores, 0)
    response=rgcca_res$call$response
    bloc_to_pred <- names(rgcca_res$call$blocks)[response]

    if (n_cores == 0)
        n_cores <- 1

    f <- quote(
        function(){

            rgcca_k <-
                set_rgcca(rgcca_res,
                          scale=scale,
                          scale_block=scale_block,
                          tol=tol,
                          scheme=scheme,
                          superblock=FALSE,
                          inds = inds,
                          method = method,
                          response=response,
                          bias = bias,
                          tau=tau,
                          ncomp=ncomp,
                          sparsity=sparsity,
                          type=type
                        ) #Rgcca on all individuals but inds
           #  
             rgcca_k_saved=rgcca_k
             rgcca_k$a <- add_variables_submodel(rgcca_res, rgcca_k$a)
             rgcca_k$astar <- add_variables_submodel(rgcca_res, rgcca_k$astar)
             rgcca_k$call$blocks <- add_variables_data(rgcca_res, rgcca_k$call$blocks)
          
             center_att <- add_variables_attr(rgcca_res, lapply(rgcca_k_saved$call$blocks, function(i) attr(i, "scaled:center")), type = "center")
          
             scale_attr <- add_variables_attr(rgcca_res, lapply(rgcca_k_saved$call$blocks, function(i) attr(i, "scaled:scale")))
            
             for (i in seq(length(rgcca_k$call$blocks))) {
                 attr(rgcca_k$call$blocks[[i]], "scaled:center") <- center_att[[i]]
                 attr(rgcca_k$call$blocks[[i]], "scaled:scale") <- scale_attr[[i]]
             }
            # Necessite les scale et les center en sortie
           respred= rgcca_predict(
                rgcca_k,
                newA = lapply(bigA, function(x) x[inds, , drop = FALSE]),
                model = model,
                fit = fit,
                bloc_to_pred = bloc_to_pred,
                # bigA = bigA,
                new_scaled = FALSE
            )

        }
    )
    if(method!="complete")
    {
        bigA <- rgcca_res$call$raw
    }
    if(method=="complete")
    {
        bigA <- intersection_list(rgcca_res$call$raw)
    }
    
  
    if (validation == "loo")
        v_inds <- seq(nrow(bigA[[1]]))
    if (validation == "kfold") {
        v_inds <- sample(nrow(bigA[[1]]))
        v_inds <- split(v_inds, sort(v_inds %% k))
    }
    if (validation == "test") {
        stop("to be implemented")
        # inds <- sample(
        #     nrow(bigA[[1]]),
        #     size = nrow(bigA[[1]]) * 0.3)
        # scores <- list(eval(f)())
        # preds <- scores$res
    }else
        {
        varlist <- c(ls(getNamespace("RGCCA")))
        # get the parameter dot-dot-dot
        args_values <- list(...)
        args_names <- names(args_values)
        n <- args_values
        if (!is.null(n))
            n <- seq(length(args_values))
        for (i in n) {
            if (!is.null(args_names[i])) {
                # dynamically asssign these values
                assign(args_names[i], args_values[[i]])
                # send them to the clusters to parallelize
                varlist <- c(varlist, args_names[i])
                # without this procedure rgcca_crossvalidation(rgcca_res, blocks = blocks2)
                # or rgcca_crossvalidation(rgcca_res, blocks = lapply(blocks, scale)
                # does not work.
            }
        }


            scores <- parallelize(
                varlist,
                seq(length(v_inds)), 
                function(i){
                    inds <- unlist(v_inds[i])
                    eval(f)()
                },
                n_cores = n_cores,
                envir = environment(),
                applyFunc = "parLapply",
                parallelization=parallelization
            ) 
        }

    list_rgcca=lapply(scores,function(x) return(x$rgcca_res))
    list_pred=lapply(scores,function(x) return(x$pred))
    list_scores=sapply(scores, function(x) x$score)
    list_res=lapply(scores, function(x) return(x$res))
    list_class.fit=lapply(scores, function(x) return(x$class.fit))
    
    if (validation %in% c("loo", "kfold")) {
        # concatenation of each test set to provide predictions for each block
        preds <- lapply(
            seq(length(rgcca_res$call$blocks)),
            function(x) Reduce(
                rbind, 
                lapply(
                    scores,
                    function(y) y$pred[[x]]
                )
            )
        )

        names(preds) <- names(rgcca_res$call$blocks)

        for (x in seq(length(preds)))
            row.names(preds[[x]]) <- row.names(bigA[[1]])
    }
     scores <- mean(unlist(lapply(scores, function(x) x$score)),na.rm=T)

    structure(
        list(scores = scores, preds = preds, rgcca_res = rgcca_res,list_scores=list_scores,list_pred=list_pred,list_rgcca=list_rgcca,list_class=list_class.fit,list_res=list_res),
        class = "cv")
}
