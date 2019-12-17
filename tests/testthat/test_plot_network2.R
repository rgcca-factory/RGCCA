#' Plot the connection between blocks (dynamic plot)
#' 
#' @inheritParams plot_ind
#' @return A dataframe with tuples of connected blocks
#' @examples
#' library(visNetwork)
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' plot_network2(rgcca_out)
#' @export
plot_network2 <- function(rgcca) {

    nodes <- get_nodes(rgcca)
    edges <- get_edges(rgcca)

    par <- ifelse("sparsity" %in% names(nodes), "sparsity", "tau")

    if (all(is.na(nodes[, par])))
        nodes[, par] <- rep("optimal", length(rgcca$blocks))

    nodes$title <- nodes$id
    nodes$label <- paste(nodes$id,
            "\nP =",
            nodes$P,
            paste0("\n", par, " ="),
            nodes[,par],
            "\nN =",
            nodes$nrow,
            sep = " ")

    edges$width <- edges$weight * 2
    nodes$color.background <- rep("#eee685", length(rgcca$blocks))

    visNetwork(
        nodes,
        edges,
        main = list(
            text = paste0("Common rows between blocks : ",
                        NROW(rgcca$blocks[[1]])),
            style = "font-family:sans;font-weight:bold;font-size:28px;text-align:center;"
        )
    ) %>%
        visNodes(
            borderWidth = 2,
            shape = "square",
            shadow = TRUE,
            color = list(
                border = "gray",
                highlight = list(background = "black", border = "darkred")
            )
        ) %>% visEdges(
            smooth = FALSE,
            shadow = TRUE,
            dashes = TRUE,
            color = list(color = "gray", highlight = "darkred")
        )

}
