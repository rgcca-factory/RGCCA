#' Cross-validation
#' 
#' Cross-validation for RGCCA
#' 
#' @inheritParams rgcca_predict
#' @inheritParams bootstrap
#' @inheritParams plot_ind
#' @param validation A character given the method for validation aong test 
#' (for test-train sets), kfold (for k-fold with k=5 by default) and loo 
#' (for leave-one-out)
#' @param k An integer given the parameter for k-fold method (5, by default)
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
rgcca_crossvalidation <- function(
    rgcca_res,
    i_block = length(rgcca_res$call$blocks),
    validation = "loo",
    type = "regression",
    fit = "lm",
    new_scaled = TRUE,
    k = 5,
    n_cores = parallel::detectCores() - 1) {

    stopifnot(is(rgcca_res, "rgcca"))
    check_blockx("i_block", i_block, rgcca_res$call$blocks)
    match.arg(validation, c("loo", "test", "kfold"))
    check_integer("k", k)
    check_integer("n_cores", n_cores, 0)

    if (n_cores == 0)
        n_cores <- 1

    if (is.null(rgcca_res$call$response))
        stop("This function requiered a RGCCA in a supervised mode.")

    bloc_to_pred <- names(rgcca_res$call$blocks)[i_block]

    f <- quote(
        function(){

            Atrain <- lapply(bigA, function(x) x[-inds, , drop = FALSE])

            if (rgcca_res$call$superblock) {
                Atrain <- Atrain[-length(Atrain)]
                rgcca_res$call$connection <- NULL
            }

            if (!is.null(rgcca_res$call$response))
                response <- length(rgcca_res$call$blocks)

            func <- quote(
                rgcca(
                    Atrain,
                    rgcca_res$call$connection,
                    response = response,
                    superblock = rgcca_res$call$superblock,
                    ncomp = rgcca_res$call$ncomp,
                    scheme = rgcca_res$call$scheme,
                    scale = FALSE,
                    type = rgcca_res$call$type,
                    verbose = FALSE,
                    init = rgcca_res$call$init,
                    bias = rgcca_res$call$bias,
                    tol = rgcca_res$call$tol,
                    method = "complete"
                )
            )

            if (rgcca_res$call$type %in% c("spls", "spca", "sgcca"))
                func$sparsity <- rgcca_res$call$c1
            else
                func$tau <- rgcca_res$call$tau

            rgcca_k <- eval(as.call(func))

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

    bigA <- intersection(rgcca_res$call$blocks)

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
        scores <- parallel::mclapply(
            seq(length(v_inds)), 
            function(i){
                inds <- unlist(v_inds[i])
                eval(f)()
            }, mc.cores = n_cores
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

    return(list(scores = scores, preds = preds))
}
