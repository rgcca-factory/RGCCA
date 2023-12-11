#' Function to postprocess the RGCCA variables
#'
#' @noRd
rgcca_postprocess_mg <- function(A, a, Y, L, g, na.rm, groups = NULL) {
  # check for parity of g
  ctrl <- all(g(-5:5) == g(5:-5))
  
  for (j in seq_along(a)) {
    if (ctrl && (a[[j]][1] < 0)) {
      a[[j]] <- -a[[j]]
      if (!is.null(groups)) {
        Y[[j]] <- pm(A[[j]], a[[j]], na.rm = na.rm)
        L[, j] <- pm(t(A[[j]]), Y[[j]], na.rm = na.rm)
      } else {
        Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
      }
    }
  }
  
  if (!is.null(groups)) {
    return(list(a = a, Y = Y, L = L))
  }
  return(list(a = a, Y = Y))
}
