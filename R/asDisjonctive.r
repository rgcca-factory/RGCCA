asDisjonctive <- function(vec, levs = NULL) {
  if (!is.null(levs)) {
    G <- factor(vec, levels = levs)
  } else {
    G <- factor(vec)
  }
  y <- data.frame(model.matrix(~ G - 1, data = G, xlev = levs))
  new_colnames <- substr(colnames(y), 2, nchar(colnames(y)))
  colnames(y) <- new_colnames
  rownames(y) <- rownames(vec)
  return(y)
}
