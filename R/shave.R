shave <- function(m, nbcomp) {
  UseMethod("shave")
}

shaving <- function(m, nbcomp) {
  UseMethod("shaving")
}

#' @export
shave.list <- function(m, nbcomp) {
  Map(shaving, m, nbcomp)
}

#' @export
shaving.matrix <- function(m, nbcomp) {
  m[, seq_len(nbcomp), drop = FALSE]
}

#' @export
shaving.double <- function(m, nbcomp) {
  m[seq_len(nbcomp)]
}

#' @export
shaving.default <- function(m, nbcomp) {
  m[seq_len(nbcomp)]
}
