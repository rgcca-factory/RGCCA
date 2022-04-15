as_disjunctive <- function(vec, levs = NULL) {
  if (NCOL(vec) > 1 || mode(vec) != "character") {
    return(vec)
  }
  if (is.null(levs) && (length(unique(vec)) == 1)) {
    stop_rgcca("Only one level in the variable to predict")
  }
  if (!is.null(levs)) {
    G <- factor(vec, levels = levs)
  } else {
    G <- factor(vec)
  }
  y <- data.frame(model.matrix(~ G - 1, data = G, xlev = levs))
  new_colnames <- substr(colnames(y), 2, nchar(colnames(y)))
  colnames(y) <- new_colnames
  rownames(y) <- rownames(vec)
  attr(y, "disjunction") <- TRUE
  return(y)
}
