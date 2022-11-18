#' Function to initialize the SGCCA variables
#'
#' @noRd
sgcca_init <- function(A, init, bias, na.rm, sparsity, response, disjunction) {
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- vapply(A, NCOL, FUN.VALUE = 1L)
  const <- sparsity * sqrt(pjs)

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

  N <- ifelse(bias, n, n - 1)
  a <- lapply(seq(J), function(j) {
    if (isTRUE(disjunction) && (j == response)) {
      a[[j]] / drop(sqrt(
        t(a[[j]]) %*% (1 / N * pm(t(A[[j]]), A[[j]], na.rm = na.rm)) %*% a[[j]]
      ))
    } else {
      soft_threshold(a[[j]], const[j])
    }
  })
  Y <- vapply(
    seq(J), function(j) pm(A[[j]], a[[j]], na.rm = na.rm),
    FUN.VALUE = double(n)
  )

  return(list(a = a, Y = Y))
}
