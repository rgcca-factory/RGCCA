#' Creates the nodes for a design matrix
#' 
#' @inheritParams plot_var_2D
#' @inheritParams select_analysis
#' @return A dataframe with rgcca$blocks in rows and the number of variables, of rows 
#' and tau or c1 in columns
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' get_nodes(rgcca = rgcca_out)
#' @export
get_nodes <- function(rgcca, tau = NULL) {

    if (is(rgcca, "sgcca")) {
        par_rgcca <- "c1"
        par.name <- "sparsity"
    } else
        par_rgcca <- par.name <- "tau"

    if (any(tau == "optimal")) {
            tau <- unlist(lapply(seq(NCOL(rgcca[[par_rgcca]])),
                function(x)
                    Reduce(paste, round(rgcca[[par_rgcca]][, x], 2))))
    }

    if (is.null(tau)) {
        if (is.matrix(rgcca[[par_rgcca]]))
            tau <-  unlist(lapply(seq(NCOL(rgcca[[par_rgcca]])),
                function(x)
                    Reduce(paste, round(rgcca[[par_rgcca]][, x], 2))))
        else
            tau <- rgcca[[par_rgcca]]
    }

    nrow <- unlist(lapply(rgcca$blocks, function(x)
            ifelse(
                is.null(attributes(x)$nrow),
                NROW(rgcca$blocks[[1]]),
                attributes(x)$nrow
            )))

    values <- list(names(rgcca$blocks), unlist(lapply(rgcca$blocks, NCOL)), nrow, tau)
    nodes <- as.data.frame(matrix(unlist(values), length(rgcca$blocks), length(values)))
    colnames(nodes) <- c("id", "P", "nrow", par.name)

    return(nodes)
}
