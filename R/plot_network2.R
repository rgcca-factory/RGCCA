#' Plot the connection between blocks (dynamic plot)
#' 
#' @inheritParams plot_ind
#' @return A dataframe with tuples of connected blocks
#' @examples
#' library(visNetwork)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks)
#' plot_network2(rgcca_out)
#' @export
plot_network2 <- function(
    rgcca_res, 
    title = paste0("Common rows between blocks : ",
                        NROW(rgcca_res$call$blocks[[1]])),
    colors =  "#eee685") {
    
    stopifnot(is(rgcca_res, "rgcca"))
    title <- paste0(title, collapse = " ")
    check_colors(colors)

    load_libraries("visNetwork")

    nodes <- get_nodes(rgcca_res)
    edges <- get_edges(rgcca_res)

    par <- ifelse("sparsity" %in% names(nodes), "sparsity", "tau")

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
    nodes$color.background <- rep(as.vector(colors[1]), length(rgcca_res$call$blocks))

    visNetwork(
        nodes,
        edges,
        main = list(
            text = title,
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
