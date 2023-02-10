#' Overlay on stats cor function to handle constant vectors.
#'
#' When a vector is constant, we say its correlation with other vectors is 0.
#' @noRd
cor2 <- function(x, y = x) {
  suppressWarnings(cor_matrix <- cor(x, y, use = "pairwise.complete.obs"))
  cor_matrix[is.na(cor_matrix)] <- 0
  return(cor_matrix)
}
