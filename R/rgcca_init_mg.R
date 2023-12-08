#' Function to initialize the RGCCA variables 
#' in the multi-group framework
#'
#' @noRd
rgcca_init_mg <- function(A, init, bias, na.rm, tau, groups = NULL) {
  J <- length(A) # number of blocks
  nb_ind <- vapply(A, NROW, FUN.VALUE = 1L) # number of individuals per group
  pjs <- vapply(A, NCOL, FUN.VALUE = 1L) # number of variables per group
  
  a <- M <- list()
  Y <- matrix(0, pjs[1], J)
  
  if (init == "svd") {
    a <- lapply(seq_along(A), function(j) {initsvd(A[[j]], dual = FALSE)})
  } else if (init == "random") {
    a <- lapply(seq_along(A), function(j) {rnorm(pjs[j])})
  }
  
  #N <- ifelse(bias, n, n - 1)
  # Only primal formulation of the problem is used
  for (j in seq_along(A)) {
    ifelse(tau[j] == 1,
           yes = {
             a[[j]] <- drop(1 / sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]]
             Y[, j] <- pm(t(A[[j]]), pm(A[[j]], a[[j]], na.rm = na.rm), na.rm = na.rm)
           },
           no = {
             M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])) * 
                              (pm(t(A[[j]]), A[[j]], na.rm = na.rm)))
             a[[j]] <- drop(1 / sqrt(t(a[[j]]) %*% M[[j]] %*% a[[j]])) *
               (M[[j]] %*% a[[j]])
             Y[, j] <- pm(t(A[[j]]), pm(A[[j]], a[[j]], na.rm = na.rm), na.rm = na.rm)
           }
    )
  }
  return(list(
    a = a, Y = Y, M = M,
    ))
}
