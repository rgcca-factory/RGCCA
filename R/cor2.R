# cor2() is similar to cor() but can handle tensor entries (flatten them). 
# @param x A numeric vector, matrix, data.frame or array.
# @param y A numeric vector, matrix, data.frame or array.
# @param forward_names A logical value. If forward_names = TRUE, and a tensor
# is flatten, its colunnames are generated from the cartesian product of its
# dimnames.
# @return \item{C}{Estimation of the correlation between x and y.}
# @title Correlation (Matrices)
# @importFrom stats cor

cor2 = function (x, y = NULL, forward_names = TRUE, ...) 
{
  if (length(dim(x)) > 2) {
    dimnames = dimnames(x)[-1]
    x        = matrix(as.vector(x), nrow = NROW(x))
    if (forward_names) {
      grid        = do.call(expand.grid, dimnames)
      colnames(x) = do.call(paste, c(grid, sep = " x "))
    }
  }
  if (is.null(y)) {
    if (is.null(dim(x))) x = as.matrix(x)
  }
  
  if (is.null(y)) {
    x = as.matrix(x)
    return(cor(x, ...))
  }
  
  if (length(dim(y)) > 2) {
    dimnames = dimnames(y)[-1]
    y        = matrix(as.vector(y), nrow = NROW(y))
    if (forward_names) {
      grid        = do.call(expand.grid, dimnames)
      colnames(y) = do.call(paste, c(grid, sep = " x "))
    }
  }
  return(cor(x, y, ...))
}
