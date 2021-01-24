#' Create either a superblock design matrix (if superblock = TRUE), or a 
#' supervised design matrix (if response != NULL) or a fully connected design 
#' matrix (if response == NULL and superblock == FALSE)
#'
#' @inheritParams rgccaNa
#' @param superblock Boolean indicating the presence of the superblock. 
#' Default = TRUE
#' @param response Position of the response block 
#' @return A binary design matrix encoding the connection between the blocks
set_connection <- function(
    blocks,
    superblock = FALSE,
    response = NULL
    ){

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
