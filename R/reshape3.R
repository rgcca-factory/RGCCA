reshape3<-function (A, n_col) 
{
  idx = rbind(seq(from = 1, to = NCOL(A) - n_col + 1, by = n_col), 
              seq(from = n_col, to = NCOL(A), by = n_col))
  B   = as.matrix(Reduce("rbind", Reduce("rbind", apply(idx, 2, function(x) list(A[, x[1]:x[2]])))))
  return(B)
}
