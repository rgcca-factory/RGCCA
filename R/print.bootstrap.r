#'Print bootstrap
#'
#'Print a bootstrap object
#'@inheritParams get_bootstrap
#'@param x A bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#'@param ... Further arguments in print
#' @return A matrix containing for each variables of each blocks, 
#' the means, 95\% intervals, bootstrap ratio, p-values and other statistics
#' \itemize{
#' \item 'estimate' for block weight vectors
#' \item 'mean' for the mean of the bootstrap block weight vectors
#' \item 'sd' for the boostrapped standard error of the block weight vectors
#' \item 'lower/upper_bound' for the lower and upper intervals from to the 
#' 'bar' parameter
#' \item 'pvals' for p-values. In the case of SGCCA, the occurrences of the 
#' bootstrap weights are distributed according to the binomial distribution 
#' while for RGCCA, their absolute value weighted by their standard deviation 
#' follows a normal distribution.
#' \item 'adjust.pval' for ajusted p-value (fdr correction by default)
#' }
#'@export
print.bootstrap = function(x, bars = "quantile", 
                           n_cores = parallel::detectCores() - 1, ...) {
    print(paste("Results based on", dim(x$bootstrap[[1]][[1]])[2], 
                "bootstrap samples"), ...)
    ncompmax = min(x$rgcca$call$ncomp)
    if(x$rgcca$call$superblock==FALSE)
    {
        for (comp in 1:ncompmax) {
            cat(paste("Dimension:", comp, "\n"))
            print(Reduce(rbind, lapply(1:length(x$rgcca$blocks), 
                                       function(block) {
                b = get_bootstrap(b = x, block = block, comp = comp, 
                                  bars = bars, display_order = FALSE, 
                                  n_cores = n_cores)
                othercols = colnames(b)[-which(colnames(b) == "estimate")]
                return(b[, c("estimate", othercols)])
            })))
        }
    }
    else
    {
        for (comp in 1:ncompmax) {
            cat(paste("Dimension:", comp, "\n"))
            print(Reduce(rbind, lapply(1:(length(x$rgcca$blocks)-1), 
                                       function(block) {
                b = get_bootstrap(b = x, block = block, comp = comp, 
                                  bars = bars, display_order = FALSE, 
                                  n_cores = n_cores)
                othercols = colnames(b)[-which(colnames(b) == "estimate")]
                return(b[, c("estimate", othercols)])
            })))
        }
    }
 
}