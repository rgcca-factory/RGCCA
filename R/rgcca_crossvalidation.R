#' Cross-validation
#' 
#' Cross-validation for RGCCA
#' 
#' @inheritParams rgcca_predict
#' @inheritParams bootstrap
#' @inheritParams plot_ind
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, response = 3)
#' rgcca_crossvalidation(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
#' rgcca_crossvalidation(rgcca_out,  validation = "test", n_cores = 1)$scores
#' rgcca_crossvalidation(rgcca_out, n_cores = 1)
#' @export
#' @seealso \link{rgcca}, \link{rgcca_predict}, \link{plot.predict}
rgcca_crossvalidation <- function(
    rgcca_res,
    response = NULL,
    validation = "loo",
    type = "regression",
    fit = "lm",
    new_scaled = TRUE,
    k = 5,
    n_cores = parallel::detectCores() - 1,
    ...) {

    stopifnot(is(rgcca_res, "rgcca"))
    if(is.null(response))
        response <- rgcca_res$call$response
    check_blockx("response", response, rgcca_res$call$blocks)
    match.arg(validation, c("loo", "test", "kfold"))
    check_integer("k", k)
    check_integer("n_cores", n_cores, 0)

    if (n_cores == 0)
        n_cores <- 1

    if (is.null(rgcca_res$call$response))
        stop("This function required an analysis in a supervised mode.")

    bloc_to_pred <- names(rgcca_res$call$blocks)[response]

    f <- quote(
        function(){

            rgcca_k <-
                set_rgcca(rgcca_res,
                          method = "complete",
                          inds = inds,
                          response = response,
                          ...)
            rgcca_k$a <- check_sign_comp(rgcca_res, rgcca_k$a)

            rgcca_predict(
                rgcca_k,
                newA = lapply(bigA, function(x) x[inds, , drop = FALSE]),
                type = type,
                fit = fit,
                bloc_to_pred = bloc_to_pred,
                bigA = bigA,
                new_scaled = TRUE
            )
        }
    )

    bigA <- attributes(
        set_rgcca(
            rgcca_res,
            method = "complete",
            inds = .Machine$integer.max,
            response = response,
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
        inds <- sample(
            nrow(bigA[[1]]),
            size = nrow(bigA[[1]]) * 0.3)
        scores <- list(eval(f)())
        preds <- scores$res
    }else{

        varlist <- c()
        # get the parameter dot-dot-dot
        args_values <- c(...)
        # get the names of the arguments of function expect the ...
        args_func_names <- names(as.list(args("rgcca_crossvalidation")))
        # get only the names of the ... args
        args_dot_names <- setdiff(names(as.list(match.call()[-1])), args_func_names)
        n <- args_values
        if(!is.null(n))
            n <- seq(length(args_values))
        for (i in n) {
            # dynamically asssign these values
            assign(args_dot_names[i], args_values[i])
            # send them to the clusters to parallelize
            varlist <- c(varlist, args_dot_names[i])
            # without this procedure rgcca_crossvalidation(rgcca_res, blocks = blocks2)
            # or rgcca_crossvalidation(rgcca_res, blocks = lapply(blocks, scale)
            # does not work.
        }

        scores <- parallelize(
            c(varlist, "check_sign_comp", "set_rgcca"),
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

    scores <- mean(unlist(lapply(scores, function(x) x$score)))

    structure(
        list(scores = scores, preds = preds, rgcca_res = rgcca_res),
        class = "cv")
}
