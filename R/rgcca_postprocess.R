rgcca_postprocess <- function(A, a, Y, g, na.rm) {
  # check for parity of g
  ctrl <- !any(g(-5:5) != g(5:-5))

  for (j in seq_along(a)) {
    if (ctrl & a[[j]][1] < 0) {
      a[[j]] <- -a[[j]]
      Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
    }
  }

  return(list(a = a, Y = Y))
}
