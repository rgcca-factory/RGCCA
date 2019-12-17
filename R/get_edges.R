#' Creates the edges for a design matrix
#' 
#' @inheritParams plot_ind
#' @return A dataframe with tuples of connected rgcca$A
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' get_edges(rgcca_out)
#' @export
get_edges <- function(rgcca) {

    J <- NCOL(rgcca$C)
    edges <- list()

    k <- 0
    for (j in seq(J)) {
        for (i in seq(J)) {
            if (i > k && rgcca$C[i, j] > 0)
                edges[[length(edges) + 1]] <-
                    c(names(rgcca$blocks)[j], names(rgcca$blocks)[i], rgcca$C[i, j])
        }
        k <- k + 1
    }

    edges <- as.data.frame(t(matrix(unlist(edges), 3, length(edges))))
    colnames(edges) <- c("from", "to", "weight")
    edges[, 3] <- as.numeric(edges[, 3])

    return(edges)
}

