# Unfold a tensor along a given mode.
# @param x  A tensor to be unfolded.
# @param mode  A mode to unfold along (default is 1).
# @return \item{x}{The resulting unfolded matrix}
# @title Tensor unfolding

unfold <- function(x, mode = 1) {
  DIM = dim(x)
  idx = c(mode, (1:length(DIM))[-mode])
  x = aperm(x, idx)
  x = matrix(as.vector(x), nrow = DIM[mode])
  return(x)
}
