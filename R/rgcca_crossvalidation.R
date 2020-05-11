#' Cross-validation
#' 
#' Uses cross-validation to validate a predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams bootstrap
#' @inheritParams plot_ind
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
    validation = "loo",
    model = "regression",
    fit = "lm",
#    new_scaled = TRUE,
    k = 5,
    n_cores = parallel::detectCores() - 1,
    ...) {

    stopifnot(is(rgcca_res, "rgcca"))
    if(is.null(rgcca_res$call$response)){
       stop("This function required an analysis in a supervised mode")}
   
    match.arg(validation, c("loo", "test", "kfold"))
    check_integer("k", k)
    check_integer("n_cores", n_cores, 0)

    if (n_cores == 0)
        n_cores <- 1

    bloc_to_pred <- names(rgcca_res$call$blocks)[length(rgcca_res$call$blocks)]

    f <- quote(
        function(){

            rgcca_k <-
                set_rgcca(rgcca_res,
                          inds = inds,
                          ...) #Rgcca on all individuals but inds
            rgcca_k$a <- check_sign_comp(rgcca_res, rgcca_k$a)

             rgcca_predict(
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

    bigA <- attributes(
        set_rgcca(
            rgcca_res,
            inds = .Machine$integer.max,
            ...
        )
    )$bigA_scaled

    if (validation == "loo")
        v_inds <- seq(nrow(bigA[[1]]))
    if (validation == "kfold") {
        v_inds <- sample(nrow(bigA[[1]]))
        v_inds <- split(v_inds, sort(v_inds %% k))
    }
    if (validation == "test") {
        print("to be implemented")
        # inds <- sample(
        #     nrow(bigA[[1]]),
        #     size = nrow(bigA[[1]]) * 0.3)
        # scores <- list(eval(f)())
        # preds <- scores$res
    }else{
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
            applyFunc = "parLapply"
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
