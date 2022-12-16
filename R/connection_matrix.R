#' Construct connection matrix based on 3 common types.
#' @noRd
connection_matrix <- function(blocks, type = "pair",
                              J = length(blocks), response = J) {
  name_blocks <- names(blocks)
  switch(type,
    "pair" = {
      connection <- 1 - diag(J)
    },
    "all" = {
      connection <- matrix(1, J, J)
    },
    "response" = {
      if (J > length(blocks)) name_blocks <- c(name_blocks, "superblock")
      connection <- matrix(0, J, J)
      connection[, response] <- connection[response, ] <- 1
      connection[response, response] <- 0
    }
  )
  rownames(connection) <- colnames(connection) <- name_blocks
  return(connection)
}
