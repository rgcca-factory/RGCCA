#' Create a matrix corresponding to a connection between the blocks
#'
#' @param blocks A list of matrix
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of
#' a superblock
#' @return A matrix corresponding to the connection between the blocks
#' @examples
#' \dontrun{
#' blocks = lapply(seq(4), function(x) matrix(runif(47 * 5), 47, 5))
#' set_connection (blocks, 'data/connection.tsv')
#' }
#' @export
set_connection <- function(
    blocks,
    superblock = FALSE
) {

    J <- length(blocks)

    if (superblock) {
        connection <- matrix(0, J, J)
        connection[seq(J - 1), J] <- connection[J, seq(J - 1)] <- 1
    } else
        connection <- 1 - diag(J)

    return(connection)
}
