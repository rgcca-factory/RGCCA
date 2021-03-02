# Helper function to generate fake data
helper.generate_blocks <- function(dims, name = FALSE) {
  J = length(dims)
  blocks = lapply(dims, function(dim) {
    block = array(rnorm(prod(dim)), dim = dim)
    if (length(dim(block)) == 1) block = as.vector(block)
    if (length(dim(block)) == 2) block = as.matrix(block)
    return(block)
  })
  if (name) {
    names(blocks) = paste0("block_", 1:J)
    for (j in 1:J) {
      if (length(dims[[j]]) == 1) {
        names(blocks[[j]]) = paste0("S_", 1:dims[[j]][1])
      } else {
        rownames(blocks[[j]]) = paste0("S_", 1:dims[[j]][1])
        if (length(dims[[j]]) == 2) {
          colnames(blocks[[j]]) = paste(
            names(blocks)[[j]], 1:dims[[j]][2], sep = "_")
        } else {
          for (d in 2:length(dims[[j]])) {
            dimnames(blocks[[j]])[[d]] = paste(
              names(blocks)[[j]], d - 1, 1:dims[[j]][d], sep = "_")
          }
        }
      }
    }
  }
  return(blocks)
}
