defl.select_tensor <- function(yy, fac, rr, nncomp, nn, nbloc)
{
  resdefl <- NULL
  pdefl <- list()
  for (q in 1:nbloc) {
    if (nn <= nncomp[q]) {
      if (length(dim(rr[[q]])) > 2) {
        W <- fac[[q]][[1]][, seq(1, nn), drop = FALSE]
        p <- diag(nrow(W)) - W %*% ginv(crossprod(W)) %*% t(W)
        p <- svd(p, nu = nrow(W) - nn)$u
        p <- cbind(p, matrix(0, nrow = nrow(W), ncol = nn))
        # p <- svd(p)$u
        resdefl[[q]] <- array(
          matrix(rr[[q]], nrow = nrow(rr[[q]])) %*%
            (diag(prod(vapply(fac[[q]][-1], nrow, FUN.VALUE = 1L))) %x% p),
          dim = dim(rr[[q]])
        )
        pdefl[[q]] <- p
      } else {
        defltmp <- deflation(rr[[q]], yy[[q]], left = FALSE)
        resdefl[[q]] <- defltmp$R
        pdefl[[q]] <- as.matrix(defltmp$p)
      }
    }
    else {
      resdefl[[q]] <- rr[[q]]
      pdefl[[q]] <- rep(0, NCOL(rr[[q]]))
    }
  }
  return(list(resdefl = resdefl, pdefl = pdefl))
}
