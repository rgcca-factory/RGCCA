#' Plot the connection between blocks
#' 
#' @inheritParams plot_ind
#' @return A dataframe with tuples of connected blocks
#' @examples
#' library(igraph)
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' plot_network(rgcca_out)
#' @export
plot_network <- function(
    rgcca, 
    title = paste0("Common rows between blocks : ",
                   NROW(rgcca$call$blocks[[1]]))) {

    stopifnot(is(rgcca, "rgcca"))
    title <- paste0(title, collapse = " ")

    # Avoid random
    set.seed(1)
    V <- E <- NULL
        
    nodes <- get_nodes(rgcca)
    edges <- get_edges(rgcca)

    par <- ifelse("sparsity" %in% names(nodes), "sparsity", "tau")

    net <- graph_from_data_frame(
        d = edges,
        vertices = nodes,
        directed = FALSE)

    if (all(is.na(nodes[, par]))) {
        nodes[, par] <- rep("optimal", length(rgcca$call$blocks))
        V(net)$call$tau <- rep(1, length(rgcca$call$blocks))
    }

    V(net)$color <- "khaki2"
    V(net)$label <- paste(
        nodes$id,
        "\nP =",
        nodes$P,
        paste0("\n", par, " ="),
        nodes[,par],
        "\nN =",
        nodes$nrow,
        sep = " ")
    V(net)$shape <- "square"
    E(net)$width <- E(net)$weight * 2
      plot(
        net,
        cex.main = 5,
        edge.color = "gray70",
        edge.lty = 2,
        vertex.frame.color = "gray50",
        vertex.label.color = "black",
        vertex.label.dist = 6,
        vertex.label.degree = 1.5,
        vertex.size = 23,
        main = title
   
    )

}
