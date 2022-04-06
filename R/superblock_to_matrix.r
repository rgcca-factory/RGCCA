superblock_to_matrix <- function(concatenedBlocks, nvar) {
  Alist <- list()
  J <- length(nvar)
  for (j in 1:J) {
    if (j == 1) {
      sel <- 1:nvar[1]
    } else {
      debut <- sum(nvar[1:(j - 1)]) + 1
      fin <- debut + (nvar[j] - 1)
      sel <- debut:fin
    }
    Alist[[j]] <- as.matrix(concatenedBlocks[, sel])
    colnames(Alist[[j]]) <- colnames(concatenedBlocks)[sel]
  }
  return(Alist)
}
