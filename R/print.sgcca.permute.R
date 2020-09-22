#' Prints the results of permutation rgcca
#' @param x A rgcca_permutation object (see  \code{\link[RGCCA]{rgcca_permutation}} )
#' @param ... Further print parameters
#' @export
#' @examples
#' data('Russett')
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' res = rgcca_permutation(A, n_run = 5, n_cores = 1)
#' print(res)
print.permutation <- function(x, ...) {

    cat("Call: ") 
    names_call <- c("type", "perm.par", "perm.value", "nperm", "quiet", "method", "tol", "scale", "scale_block", "superblock")
    char_to_print <- ""
    for (name in names_call) {

        if (name == "ncomp") {
            if (length(x$call$ncomp) > 1) {
                value <- (paste(x$call$ncomp, sep = "", collapse = ","))
                value <- paste0("c(", value, ")")
            }
        }
        if (name != "ncomp") {
            value <- x$call[[name]]
        }

        quo <- ifelse(is.character(value) & name != "ncomp", "'", "")
        vir <- ifelse(name == names_call[length(names_call)], "", ", ")
        char_to_print <- paste(char_to_print, name, "=", quo, value, quo, vir, collapse = "", sep = "")
    }
    cat(char_to_print, "\n")

    cat("The design matrix is:\n")
    colnames(x$call$connection) <- rownames(x$call$connection) <- names(x$a)
    print(x$call$connection)
    cat("\n")

    if (is.function(x$call$scheme)) {
        cat("The", deparse(x$call$scheme), "scheme was used.", fill = TRUE)
    } else {
        cat("The", x$call$scheme, "scheme was used.", fill = TRUE)
    }

    c1s <- round(x$penalties, 4)
    rownames(c1s) <- 1:NROW(c1s)
    cat(fill = TRUE)
    cat("Tuning parameters used: ", fill = TRUE)
    print(c1s, quote = FALSE, ...)
    cat("\n")

    tab <- round(cbind(x$crit, x$means, x$sds, x$zstat, x$pvals), 3)
    dimnames(tab) <- list(paste("Tuning parameter set ", sep = "", 1:length(x$pvals)), c("Crit", "Crit Perm", "Sd", "Z", "P-Value"))
    print(tab, quote = FALSE, ...)

    cat("Tuning parameters corresponding to highest z score: \n")
    cat(paste(round(x$bestpenalties, 3), collapse = ", "), "\n")
    cat("Highest z score: ", max(x$zstat), "\n")
    cat("P-value corresponding to highest z score: ", x$pvals[which.max(x$zstat)], fill = TRUE)
}