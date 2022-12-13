#' Function to postprocess the SGCCA variables
#'
#' @noRd
sgcca_postprocess <- function(
    A, a, Y, g, na.rm, sparsity, tol, response, disjunction
  ) {
  pjs <- vapply(A, NCOL, FUN.VALUE = 1L)

  # check for parity of g
  ctrl <- all(g(-5:5) == g(5:-5))

  for (j in seq_along(a)) {
    if (ctrl && (a[[j]][1] < 0)) {
      a[[j]] <- -a[[j]]
      Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
    }
  }

  l2_sat <- vapply(a, function(x) norm(x, "2"), FUN.VALUE = 1.0)
  if (disjunction) {
    l2_sat <- l2_sat[-response]
  }
  if (max(abs(l2_sat - 1)) > tol) {
    for (i in which(abs(l2_sat - 1) > tol)) {
      if (l2_sat[i] < .Machine$double.eps) {
        warning(
          "Norm2 of the block weight vector #",
          i, " is too small :", l2_sat[i]
        )
      } else {
        nMAX <- length(which(a[[i]] != 0))
        warning(
          "The l2 constraint is not saturated for block #", i,
          ". The sparsity parameter has to be in the range [",
          sqrt(nMAX / pjs[i]),
          ", 1] and is equal to ", sparsity[i], "."
        )
      }
    }
  }

  return(list(a = a, Y = Y))
}
