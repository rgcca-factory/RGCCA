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
#' rgcca_out = rgcca.analyze(blocks)
#' rgcca_crossvalidation(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
#' rgcca_crossvalidation(rgcca_out,  validation = "test", n_cores = 1)$scores
#' rgcca_crossvalidation(rgcca_out, n_cores = 1)
#' @export
rgcca_crossvalidation <- function(
    rgcca,
    i_block = length(rgcca$blocks),
    validation = "loo",
    type = "regression",
    fit = "lm",
    new_scaled = TRUE,
    k = 5,
    n_cores = parallel::detectCores() - 1) {
    
    bloc_to_pred = names(rgcca$blocks)[i_block]

    match.arg(validation, c("test", "kfold", "loo"))

    f <- quote(
        function(){

            Atrain <- lapply(bigA, function(x) x[-inds, ])

            if (is(rgcca, "sgcca"))
                tau <- rgcca$c1
            else
                tau <- rgcca$tau

            if (rgcca$superblock) {
                Atrain <- Atrain[-length(Atrain)]
                rgcca$C <- NULL
            }

            rgcca_k <- rgcca.analyze(
                Atrain,
                rgcca$C,
                superblock = rgcca$superblock,
                tau = tau,
                ncomp = rgcca$ncomp,
                scheme = rgcca$scheme,
                scale = FALSE,
                type = class(rgcca),
                verbose = FALSE,
                init = rgcca$init,
                bias = rgcca$bias,
                tol = rgcca$tol,
                method="complete"
            )

            rgcca_k$a <- check_sign_comp(rgcca, rgcca_k$a)

            rgcca_predict(
                rgcca_k,
                newA = lapply(bigA, function(x) x[inds, ]),
                type = type,
                fit = fit,
                bloc_to_pred = bloc_to_pred,
                bigA = bigA,
                new_scaled = TRUE
            )
        }
    )

    bigA <- rgcca$blocks

    if (validation == "loo")
        v_inds <- seq(nrow(rgcca$blocks[[1]]))
    if (validation == "kfold") {
        v_inds <- sample(nrow(rgcca$blocks[[1]]))
        v_inds <- split(v_inds, sort(v_inds %% k))
    }

    if (validation == "test") {
        inds <- sample(
            nrow(rgcca$blocks[[1]]),
            size = nrow(rgcca$blocks[[1]]) * 0.3)
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
            seq(length(rgcca$blocks)),
            function(x) Reduce(
                rbind, 
                lapply(
                    scores,
                    function(y) y$pred[[x]]
                    )
                )
            )

        names(preds) <- names(rgcca$blocks)

    for (x in seq(length(preds)))
        row.names(preds[[x]]) <- row.names(rgcca$blocks[[1]])

    }

    scores <- mean(unlist(lapply(scores, function(x) x$score)))

    return(list(scores = scores, preds = preds))
}
