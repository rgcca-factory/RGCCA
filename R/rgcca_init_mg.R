#' Function to initialize the RGCCA variables 
#' in the multi-group framework
#'
#' @noRd
rgcca_init_mg <- function(A, init, na.rm, tau, groups = NULL) {
  J <- length(A) # number of groups
  nb_ind <- vapply(A, NROW, FUN.VALUE = 1L) # number of individuals per group
  p <- NCOL(A[[1]]) # number of variables per group
  
  a <- M <- Y <- list()
  L <- matrix(0, p, J)
  
  if (bias) {
    N <- nb_ind
  } else {
    N <- nb_ind - 1
  }
  
  if (init == "svd") {
    a <- lapply(seq_along(A), function(j) {initsvd(A[[j]], dual = FALSE)})
  } else if (init == "random") {
    a <- lapply(seq_along(A), function(j) {rnorm(p)})
  }
  
  # Only primal formulation of the problem is used
  for (j in seq_along(A)) {
    ifelse(tau[j] == 1,
           yes = {
             a[[j]] <- drop(1 / sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]]
             Y[[j]] <- pm(A[[j]], a[[j]], na.rm = na.rm)
             L[, j] <- pm(t(A[[j]]), Y[[j]], na.rm = na.rm)
           },
           no = {
             M[[j]] <- ginv(tau[j] * diag(p) + ((1 - tau[j])) * 1 / N[[j]] *
                              (pm(t(A[[j]]), A[[j]], na.rm = na.rm)))
             a[[j]] <- drop(1 / sqrt(t(a[[j]]) %*% M[[j]] %*% a[[j]])) *
               (M[[j]] %*% a[[j]])
             Y[[j]] <- pm(A[[j]], a[[j]], na.rm = na.rm)
             L[, j] <- pm(t(A[[j]]), Y[[j]], na.rm = na.rm)
           }
    )
  }
  return(list(
    a = a, Y = Y, L = L, M = M
  ))
}
