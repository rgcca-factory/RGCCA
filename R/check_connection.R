#' Check the format of the connection matrix
#'
#' @param c A symmetric matrix containing 1 and 0
#' @inheritParams set_connection
check_connection <- function(c, blocks) {

    msg <- "The connection file should"

    if (!isSymmetric.matrix(unname(c)))
        stop(paste(msg, "be a symmetric matrix."), exit_code = 103)

    # d <- unique(diag(c))
    # if (length(d) != 1 || d != 0)
    #     stop("The diagonal of the connection matrix file should be 0.",
    #         exit_code = 105)

    x <- unique(c %in% c(0, 1))
    if (length(x) != 1 || x != TRUE)
        stop(paste(msg, "contain only 0 or 1."), exit_code = 106)

    if (all(c == 0))
        stop(paste(msg, "not contain only 0."), exit_code = 107)

    invisible(check_size_blocks(blocks, "connection matrix", c))

    # TODO: warning if superblock = TRUE
}
