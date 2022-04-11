#' Function to initialize the SGCCA variables
#'
#' @noRd
sgcca_init <- function(A, init, bias, na.rm, const, pjs, J, n) {
  a <- list()
  Y <- matrix(0, n, J)

  #  Choose J arbitrary vectors
  if (init == "svd") {
    # SVD initialization for a_j
    a <- lapply(A, function(x) {
      return(initsvd(x, dual = FALSE))
    })
    a <- lapply(a, function(x) {
      return(as.vector(x))
    })
  } else if (init == "random") {
    a <- lapply(pjs, rnorm)
  }

  a <- lapply(seq(J), function(b) soft_threshold(a[[b]], const[b]))
  Y <- vapply(
    seq(J), function(b) pm(A[[b]], a[[b]], na.rm = na.rm),
    FUN.VALUE = double(n)
  )

  return(list(a = a, Y = Y))
}
