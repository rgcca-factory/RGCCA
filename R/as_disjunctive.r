#' One hot encoding of a vector
#'
#' This function tries to one hot encode a vector of characters or leave it
#' untouched if not relevant.
#' @param vec vector to transform
#' @param levs factor levels associated to the possible values of \code{vec}.
#' If NULL, levels are directly taken from the values of \code{vec}.
#' @return One hot encoded version of \code{vec} or \code{vec} unchanged.
#' @noRd
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
  # Change na_option locally to keep rows of NA
  current_na_action <- options("na.action")
  options(na.action = "na.pass")
  y <- data.frame(model.matrix(~ G - 1, data = G, xlev = levs))
  options(current_na_action)

  new_colnames <- substr(colnames(y), 2, nchar(colnames(y)))
  colnames(y) <- new_colnames
  rownames(y) <- rownames(vec)
  y <- y[, -1, drop = FALSE]
  return(y)
}
