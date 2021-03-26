#'Print bootstrap
#'
#'Print a bootstrap object
#' @param x A fitted bootstrap object (see \code{\link[RGCCA]{bootstrap}})
#' @param type Character string indicating the bootstrapped object to print:
#' block-weight vectors ("weight", default) or block-loading vectors
#' ("loadings").
#'@param ... Further arguments in print
#' @return A matrix containing for each variables of each blocks,
#' the means, 95\% intervals, bootstrap ratio, p-values and other statistics.
#' \itemize{
#' \item 'estimate' for block weight/loading vectors.
#' \item 'mean' for the mean of the bootstrap block weight/loading vectors.
#' \item 'sd' for the boostrapped estimate of the standard deviation of the
#' block weight/loading vectors.
#' \item 'lower/upper_bound' for the lower and upper intervals.
#' \item 'pvals' for p-values. In the case of SGCCA, the occurrences of the
#' bootstrap weights are distributed according to the binomial distribution
#' while for RGCCA, their absolute value divided by their standard deviation
#' follows a normal distribution.
#' \item 'adjust.pval' for ajusted p-value (fdr correction by default)
#' }
#'@export
print.bootstrap = function(x, type = "weight", ...) {
    print(paste0("Extract statistics on the block-", type, " vectors from ",
                NCOL(x$bootstrap[[1]][[1]][[1]]), " bootstrap samples"), ...)

    ncompmax = min(x$rgcca$call$ncomp)

    if(x$rgcca$call$superblock==FALSE)
    {
        for (comp in 1:ncompmax) {
            cat(paste("Dimension:", comp, "\n"))
            print(Reduce(rbind, lapply(1:length(x$rgcca$call$blocks),
                                       function(block) {
                b = get_bootstrap(b = x, type = type, block = block,
                                  comp = comp, display_order = FALSE)
                othercols = colnames(b)[-which(colnames(b) == "estimate")]
                return(b[, c("estimate", othercols)])
            })))
        }
    }
    else
    {
        for (comp in 1:ncompmax) {
            cat(paste("Dimension:", comp, "\n"))
            print(Reduce(rbind, lapply(1:(length(x$rgcca$call$blocks)-1),
                                       function(block) {
                b = get_bootstrap(b = x, type = type, block = block,
                                  comp = comp, display_order = FALSE)
                othercols = colnames(b)[-which(colnames(b) == "estimate")]
                return(b[, c("estimate", othercols)])
            })))
        }
    }
}
