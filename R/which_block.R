which_block <- function(x) {
  l <- vector("list", max(x))
  l[[1]] <- seq_along(x)
  if (max(x) > 1) {
    for (i in seq(2, max(x))) {
      l[[i]] <- which(x - (i - 1) * rep(1, length(x)) > 0)
    }
  }
  return(l)
}
