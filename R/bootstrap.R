#' Compute bootstrap
#'
#' Computing boostrap of RGCCA
#'
#' @inheritParams rgcca
#' @inheritParams plot_var_2D
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
bootstrap <- function(
    rgcca,
    n_boot = 5,
    scale = FALSE,
    n_cores = parallel::detectCores() - 1,
    ...) {

    # TODO : nboot > 1
    stopifnot(!missing(rgcca))

    if (n_cores == 0)
        n_cores <- 1

    # if (any(unlist(lapply(rgcca$call$blocks, NCOL) > 1000)))
    #     verbose <- TRUE

    cat("Bootstrap in progress...")

    W <- parallel::mclapply(
        seq(n_boot), 
        function(x) bootstrap_k(rgcca, ...), 
        mc.cores = n_cores)

    cat("OK.\n", append = TRUE)
    
  
    return(list(W=W,call=rgcca))
}
