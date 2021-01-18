#' @title Print the call of rgcca_permutation results
#' @param x A permutation object (see \code{\link{rgcca_permutation}})
#' @param ... other parameters used in print (for the displaying of matrices)
#' @examples
#' data("Russett")
#' A <- list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' perm <- rgcca_permutation(A, par_type="tau",n_run = 2, n_cores = 1)
#' print(perm)
#'@export
print.permutation <- function(x,...) {
    tab_res <- cbind(x$crit, x$means, x$sds, x$zstat, x$pvals)
    colnames(tab_res) <- c("RGCCA crit", "Perm. crit", "S.D.", "Z", "P-value")
    print(round(cbind(x$penalties, tab_res), 3),...)
}
