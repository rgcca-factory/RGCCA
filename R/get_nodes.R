# Creates the nodes for a design matrix
#
# @inheritParams plot_var_2D
# @return A dataframe with rgcca_res$call$blocks in rows and the number of variables, of rows
# and tau or sparsity in columns
get_nodes <- function(rgcca_res) {

    if ( tolower(rgcca_res$call$type) %in% c("sgcca", "spls", "spca")) {
        par_rgcca <- "sparsity"
        par_name <- "sparsity"
    } else
        par_rgcca <- par_name <- "tau"

    if (is.matrix(rgcca_res$call[[par_rgcca]]))
        penalty <-  unlist(
            lapply(
                seq(NCOL(rgcca_res$call[[par_rgcca]])),
            function(x)
                Reduce(paste, round(rgcca_res$call[[par_rgcca]][, x], 2))))
    else
        penalty <- round(rgcca_res$call[[par_rgcca]], 2)

    nrow <- unlist(lapply(rgcca_res$call$blocks, function(x)
            ifelse(
                is.null(attributes(x)$nrow),
                NROW(rgcca_res$call$blocks[[1]]),
                attributes(x)$nrow
            )))

    values <- list(names(rgcca_res$call$blocks), unlist(lapply(rgcca_res$call$blocks, NCOL)), nrow, penalty)
    nodes <- as.data.frame(matrix(unlist(values), length(rgcca_res$call$blocks), length(values)))
    colnames(nodes) <- c("id", "P", "nrow", par_name)

    return(nodes)
}
