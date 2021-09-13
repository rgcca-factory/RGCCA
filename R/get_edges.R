# Creates the edges for a design matrix
#
# @inheritParams plot_ind
# @return A dataframe with tuples of connected rgcca_res$call$blocks



get_edges <- function(rgcca_res) {

    J <- NCOL(rgcca_res$call$connection)
    edges <- list()

    k <- 0
    for (j in seq(J)) {
        for (i in seq(J)) {
            if (i > k && rgcca_res$call$connection[i, j] > 0)
                edges[[length(edges) + 1]] <-
                    c(names(rgcca_res$call$blocks)[j], names(rgcca_res$call$blocks)[i], rgcca_res$call$connection[i, j])
        }
        k <- k + 1
    }

    edges <- as.data.frame(t(matrix(unlist(edges), 3, length(edges))))
    colnames(edges) <- c("from", "to", "weight")
    edges[, 3] <- as.numeric(edges[, 3])

    return(edges)
}
