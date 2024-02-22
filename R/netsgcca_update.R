#' Function to perform netSGCCA update
#'
#' @noRd
netsgcca_update <- function(A, bias, na.rm, sparsity, lambda,
                         # group_sparsity, 
                         graph_laplacians, 
                         response, disjunction,
                         dg, C, a, Y, init_object) {
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- vapply(A, NCOL, FUN.VALUE = 1L)
  const <- sparsity * sqrt(pjs)

  Z <- matrix(0, n, J)
  N <- ifelse(bias, n, n - 1)
  penalty_blocs <- which(lambda != 0)
  
  for (j in seq_len(J)) {
    dgx <- dg(cov2(Y[, j], Y, bias = bias))
    CbyCovq <- drop(C[j, ] * dgx)
    Z[, j] <- Y %*% CbyCovq
    
    if (j %in% penalty_blocs) {
      grad <- pm(t(A[[j]]), Z[, j], na.rm = na.rm) + 
        lambda[j] * graph_laplacians[[j]] %*% a[[j]]
    } else {
      grad <- pm(t(A[[j]]), Z[, j], na.rm = na.rm)
    }
    
    if (all(grad == 0)) {
      a[[j]] <- a[[j]] * 0
    } else {
      if (disjunction && (j == response)) {
        M <- ginv(1 / N * pm(t(A[[j]]), A[[j]], na.rm = na.rm))
        a[[j]] <- M %*% grad / drop(sqrt(t(grad) %*% M %*% grad))
      } else {
        a[[j]] <- soft_threshold(grad, const[j])
      }
    }
    Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
  }

  return(list(a = a, Y = Y))
}
