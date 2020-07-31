#' Create a matrix corresponding to a connection between the blocks
#'
#' @param blocks A list of matrix
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of
#' a superblock
#' @param response if not NULL, number corresponding to the response block
#' @return A matrix corresponding to the connection between the blocks


set_connection <- function(
    blocks,
    superblock = FALSE,response=NULL
) {

    J <- length(blocks)

    if (superblock) {
        connection <- matrix(0, J, J)
        connection[seq(J - 1), J] <- connection[J, seq(J - 1)] <- 1
    }
    else{
        if(!is.null(response))
        {
            connection <- matrix(0, J, J)
            Resp=response
            notResp=(1:J)[-Resp]
            connection[notResp, Resp] <- connection[Resp, notResp] <- 1
        } 
        else
        {
            connection <- 1 - diag(J)
        }
            
    }
   
    row.names(connection) <- names(blocks) -> colnames(connection)
    return(connection)
}
