#' Compute bootstrap
#'
#' Computing boostrap of RGCCA in order to visualize the stability of the weights found in RGCCA
#' @param rgcca_res Result of a RGCCA (see  \code{\link[RGCCA]{rgcca}} )
#' @param n_boot A integer for the number of boostrap
#' @param n_cores An integer for the number of cores used in parallelization 
#' @param ... other RGCCA parameters # TODO
#' @return A list of RGCCA bootstrap weights
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
#' bootstrap(rgcca_out, n_boot = 2, n_cores = 1, blocks = lapply(blocks, scale),
#'  superblock = FALSE)
#' @export
#' @seealso \code{\link[RGCCA]{plot.bootstrap}} , \code{\link[RGCCA]{print.bootstrap}} 
bootstrap <- function(
    rgcca_res,
    n_boot = 5,
    n_cores = parallel::detectCores() - 1,
    ...) {

    stopifnot(is(rgcca_res, "rgcca"))
    check_integer("n_boot", n_boot)
    check_integer("n_cores", n_cores, 0)

    if (n_cores == 0)
        n_cores <- 1

    # if (any(unlist(lapply(rgcca$call$blocks, NCOL) > 1000)))
    #     verbose <- TRUE

    cat("Bootstrap in progress...")

    blocks <- NULL

    varlist <- c()
    # get the parameter dot-dot-dot
    args_values <- c(...)
    # get the names of the arguments of function expect the ...
    args_func_names <- names(as.list(args("bootstrap")))
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

    W <- parallelize(
        c(varlist, "bootstrap_k", "remove_null_sd", "check_sign_comp", "set_rgcca", "blocks"),
        seq(n_boot), 
        function(x) bootstrap_k(rgcca_res, ...), 
        n_cores = n_cores,
        envir = environment(),
        applyFunc = "parLapply")

    cat("OK.\n", append = TRUE)

    return(structure(list(bootstrap = W, rgcca = rgcca_res), class = "bootstrap"))
}
