#' Creates the nodes for a design matrix
#' 
#' @inheritParams plot_var_2D
#' @inheritParams select_analysis
#' @param penalty NULL, "optimal", "sparsity" or "tau"
#' @return A dataframe with rgcca$call$blocks in rows and the number of variables, of rows
#' and tau or sparsity in columns
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' get_nodes(rgcca = rgcca_out)
#' @export
get_nodes <- function(rgcca, penalty = NULL) {

    if ( rgcca$call$type %in% c("sgcca", "spls", "spca")) {
        par_rgcca <- "sparsity"
        par.name <- "sparsity"
    } else
        par_rgcca <- par.name <- "tau"

    if (any(penalty == "optimal")) {
            penalty <- unlist(lapply(seq(NCOL(rgcca$call[[par_rgcca]])),
                function(x)
                    Reduce(paste, round(rgcca$call[[par_rgcca]][, x], 2))))
    }

    if (is.null(penalty)) {
        if (is.matrix(rgcca$call[[par_rgcca]]))
            penalty<-  unlist(lapply(seq(NCOL(rgcca$call[[par_rgcca]])),
                function(x)
                    Reduce(paste, round(rgcca$call[[par_rgcca]][, x], 2))))
        else
            penalty<- rgcca$call[[par_rgcca]]
    }

    nrow <- unlist(lapply(rgcca$call$blocks, function(x)
            ifelse(
                is.null(attributes(x)$nrow),
                NROW(rgcca$call$blocks[[1]]),
                attributes(x)$nrow
            )))

    values <- list(names(rgcca$call$blocks), unlist(lapply(rgcca$call$blocks, NCOL)), nrow, penalty)
    nodes <- as.data.frame(matrix(unlist(values), length(rgcca$call$blocks), length(values)))
    colnames(nodes) <- c("id", "P", "nrow", par.name)

    return(nodes)
}
