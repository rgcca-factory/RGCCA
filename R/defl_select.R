#' The function defl_select() computes residual matrices \eqn{X_{1,h+1},
#' \ldots, X_{J,h+1}}. These residual matrices are determined according to the
#' following formula: \eqn{X_{j,h+1} = X_{jh} - y_{jh} p_{jh}^t}.
#' @inheritParams rgccad
#' @param yy  A matrix that contains the GCCA block components of each block:
#' \eqn{y_{1h}, \ldots, y_{Jh}}
#' @param rr  A list that contains the residual matrices
#' \eqn{X_{1h}, \ldots, X_{Jh}}
#' @param nncomp A \eqn{1 \times J} vector that contains the number of
#' components to compute for each block.
#' @param nn  A \eqn{1 \times J} vector that contains the numbers of
#' deflations for each block
#' @param nbloc Number of blocks.
#' @return \item{resdefl}{A list of \eqn{J} elements that contains
#' \eqn{X_{1,h+1}, \ldots, X_{J,h+1}}.}
#' @return \item{pdefl}{A list of \eqn{J} elements that contains
#' \eqn{p_{1h}, \ldots, p_{Jh}}.}
#' @title deflation function
#' @noRd
defl_select <- function(yy, rr, nncomp, nn, nbloc, na.rm = TRUE,
                        response = NULL, left = TRUE) {
  resdefl <- NULL
  pdefl <- list()
  for (q in seq(nbloc)) {
    is_response <- !is.null(response) && (q == response)
    if ((nn <= nncomp[q]) && !is_response) {
      defltmp <- deflation(as.matrix(rr[[q]]),
        as.matrix(yy[[q]]),
        na.rm = na.rm, left = left
      )
      resdefl[[q]] <- defltmp$R
      pdefl[[q]] <- as.matrix(defltmp$p)
    } else {
      resdefl[[q]] <- rr[[q]]
      pdefl[[q]] <- rep(0, NCOL(rr[[q]]))
    }
  }
  return(list(resdefl = resdefl, pdefl = pdefl))
}
