#' Check the format of the connection matrix
#'
#' @inheritParams rgccad
#' @inheritParams set_connection
check_connection <- function(C, blocks) {

    msg <- "The connection file should"

    if (!isSymmetric.matrix(unname(C)))
        stop_rgcca(paste(msg, "be a symmetric matrix."), exit_code = 103)

    # d <- unique(diag(C))
    # if (length(d) != 1 || d != 0)
    #     stop_rgcca("The diagonal of the connection matrix file should be 0.",
    #         exit_code = 105)

    x <- unique(C %in% c(0, 1))
    if (length(x) != 1 || x != TRUE)
        stop_rgcca(paste(msg, "contain only 0 or 1."), exit_code = 106)

    if (all(C == 0))
        stop_rgcca(paste(msg, "not contain only 0."), exit_code = 107)

    if(is.null(rownames(C)) || is.null(colnames(C)))
        rownames(C) <- names(blocks) -> colnames(C)

    if (!all(rownames(C) %in% names(blocks)) || 
        !all(colnames(C) %in% names(blocks)))
        stop_rgcca(paste(msg,
            "have the rownames and the colnames that match with the names of the blocks."),
            exit_code = 108)

    invisible(check_size_blocks(blocks, "connection matrix", C))
    
    return(C)
    # TODO: warning if superblock = TRUE
}
