#' Cross-validation
#' 
#' Uses cross-validation to validate a predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams rgcca
#' @inheritParams bootstrap
#' @inheritParams plot_ind
#' @param k An integer giving the number of folds (if validation = 'kfold').
#' @param validation A character for the type of validation among "loo", "kfold"
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks, response = 3,superblock=FALSE)
#' res=rgcca_crossvalidation(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
#' rgcca_crossvalidation(rgcca_out, n_cores = 1)
#' rgcca_crossvalidation(rgcca_out, validation = "loo", k = 5, fit = "r2", n_cores = 1)$scores 
#' # 0.6282727
#' @export
#' @seealso \link{rgcca}, \link{rgcca_predict}, \link{plot.predict}
rgcca_crossvalidation <- function(
    rgcca_res,
    validation = "kfold",
    model = "regression",
    fit = "lm",
    regress_on = "block",
    new_scaled = TRUE,
    k = 5,
    blocks = NULL,
    scale = NULL,
    scale_block = NULL,
    tol = 1e-8,
    scheme = NULL,
    method = NULL,
    type = NULL,
    init = NULL,
    bias = NULL,
    connection = NULL,
    ncomp = NULL,
    response = NULL,
    tau = NULL,
    sparsity = NULL,
    n_cores = parallel::detectCores() - 1,
    parallelization = TRUE, 
    ...) {

    # if (!is.null(rgcca_res)) {
        stopifnot(is(rgcca_res, "rgcca"))
        # message("All parameters were imported from a rgcca object.")
        blocks <- rgcca_res$call$blocks
        scale_block <- rgcca_res$call$scale_block
        scale <- rgcca_res$call$scale
        scheme <- rgcca_res$call$scheme
        response <- rgcca_res$call$response
        tol <- rgcca_res$call$tol
        method <- rgcca_res$call$method
        init <- rgcca_res$call$init
        bias <- rgcca_res$call$bias
        blocks <- rgcca_res$call$raw
        connection <- rgcca_res$call$connection
        tau <- rgcca_res$call$tau
        ncomp <- rgcca_res$call$ncomp
        sparsity <- rgcca_res$call$sparsity
    # }

    if (is.null(response))
        stop_rgcca("This function required an analysis in a supervised mode")
    if (!is.null(parallelization))
        check_boolean("parallelization", parallelization)
    match.arg(validation, c("loo", "kfold"))
    check_integer("k", k, min = 2)
    check_integer("n_cores", n_cores, min = 0)
    bloc_to_pred <- names(blocks)[response]

    if (n_cores == 0)
        n_cores <- 1

    f <- quote(
        function(){

            rgcca_k <- set_rgcca(
                rgcca_res,
                scale = scale,
                scale_block = scale_block,
                tol = tol,
                scheme = scheme,
                superblock = FALSE,
                inds = inds,
                method = method,
                response = response,
                bias = bias,
                tau = tau,
                ncomp = ncomp,
                sparsity = sparsity,
                type = type
            )
            rgcca_k_saved <- rgcca_k
            rgcca_k$a <- add_variables_submodel(rgcca_res, rgcca_k$a)
            rgcca_k$astar <- add_variables_submodel(rgcca_res, rgcca_k$astar)
            rgcca_k$call$blocks <- add_variables_data(rgcca_res, rgcca_k$call$blocks)

            center_att <- add_variables_attr(
                rgcca_res, 
                lapply(
                    rgcca_k_saved$call$blocks, 
                    function(i) attr(i, "scaled:center")), 
                type = "center")
            scale_attr <- add_variables_attr(
                rgcca_res, 
                lapply(
                    rgcca_k_saved$call$blocks, 
                    function(i) attr(i, "scaled:scale")))

            for (i in seq(length(rgcca_k$call$blocks))) {
                attr(rgcca_k$call$blocks[[i]], "scaled:center") <- center_att[[i]]
                attr(rgcca_k$call$blocks[[i]], "scaled:scale") <- scale_attr[[i]]
            }
            respred <- rgcca_predict(
                rgcca_k,
                newA = lapply(bigA, function(x) x[inds, , drop = FALSE]),
                model = model,
                fit = fit,
                bloc_to_pred = bloc_to_pred,
                new_scaled = FALSE,
                regress_on = regress_on
            )
        }
    )

    if (method != "complete") 
        bigA <- rgcca_res$call$raw 
    else
        bigA <- intersection_list(rgcca_res$call$raw)
  
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
    }else {

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
            parallelization = parallelization
        ) 
    }

    if (validation %in% c("loo", "kfold")) {
        # concatenation of each test set to provide predictions for each block
        preds <- lapply(
            seq(length(blocks)),
            function(x) Reduce(
                rbind, 
                lapply(
                    scores,
                    function(y) y$pred[[x]]
                )
            )
        )

        names(preds) <- names(blocks)

        for (x in seq(length(preds)))
            row.names(preds[[x]]) <- row.names(bigA[[1]])
    }

    structure(
        list(
            scores = mean(unlist(lapply(scores, function(x) x$score)), na.rm = T),
            preds = preds,
            rgcca_res = rgcca_res,
            list_scores = sapply(scores, function(x) x$score),
            list_pred = lapply(scores, function(x) return(x$pred)),
            list_rgcca = lapply(scores, function(x) return(x$rgcca_res)),
            list_class = lapply(scores, function(x) return(x$class.fit)),
            list_res = lapply(scores, function(x) return(x$res))
        ),
        class = "cv"
    )
}
