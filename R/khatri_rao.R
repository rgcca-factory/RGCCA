# khatri_rao computes the Khatri-Rao products between two matrices
# @param x A numeric matrix.
# @param y A numeric matrix.
# @return \item{res}{Khatri-Rao product between x and y.}
# @title Khatri-Rao product

khatri_rao <- function(x, y){
  if ((length(dim(x)) != 2) || (length(dim(y)) != 2)){
    stop_rgcca("x and y must be matrices")
  }
  ncol   = dim(x)[2]
  if (ncol != dim(y)[2]){
    stop_rgcca("x and y must have the same number of columns")
  }
  res = sapply(1:ncol, function(i) x[, i] %x% y[, i])
  return(res)
}