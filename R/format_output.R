#' Format output of rgcca function
#'
#' @noRd
format_output <- function(func_out, rgcca_args, opt, blocks, disjunction) {
  ### Compute AVE
  blocks_AVE <- seq_along(blocks)
  names_AVE <- names(blocks)
  ncomp_AVE <- rgcca_args$ncomp
  pjs <- vapply(blocks, NCOL, FUN.VALUE = 1L)
  if (isTRUE(disjunction)) {
    blocks_AVE <- blocks_AVE[-rgcca_args$response]
    names_AVE <- names_AVE[-rgcca_args$response]
    ncomp_AVE <- ncomp_AVE[-rgcca_args$response]
    pjs <- pjs[-rgcca_args$response]
  }

  AVE_inner <- vapply(seq(max(rgcca_args$ncomp)), function(n) {
    sum(rgcca_args$connection * cor(
      do.call(cbind, lapply(func_out$Y, function(y) y[, n]))
    )^2 / 2) / (sum(rgcca_args$connection) / 2)
  }, FUN.VALUE = double(1L))

  AVE_X <- lapply(blocks_AVE, function(j) {
    apply(func_out$Y[[j]], 2, rsq, blocks[[j]])
  })
  AVE_X_cum <- lapply(blocks_AVE, function(j) {
    vapply(
      seq_len(NCOL(func_out$Y[[j]])),
      function(p) rsq(func_out$Y[[j]][, seq(p)], blocks[[j]]),
      FUN.VALUE = 1.0
    )
  })
  AVE_X_cor <- lapply(AVE_X_cum, function(x) {
    y <- c(0, x[-length(x)])
    x - y
  })

  outer <- matrix(unlist(AVE_X_cor), nrow = max(ncomp_AVE))
  AVE_outer <- as.vector((outer %*% pjs) / sum(pjs))
  AVE_X <- shave(AVE_X, ncomp_AVE)
  AVE_X_cor <- shave(AVE_X_cor, ncomp_AVE)
  func_out$AVE <- list(
    AVE_X = AVE_X, AVE_X_cor = AVE_X_cor,
    AVE_outer = AVE_outer, AVE_inner = AVE_inner
  )
  func_out$AVE_inner <- NULL

  names(func_out$AVE$AVE_X) <- names_AVE
  names(func_out$AVE$AVE_X_cor) <- names_AVE

  ### Set names and shave
  for (j in seq_along(blocks)) {
    rownames(func_out$a[[j]]) <- colnames(blocks[[j]])
    rownames(func_out$Y[[j]]) <- rownames(blocks[[j]])
    colnames(func_out$Y[[j]]) <- paste0("comp", seq_len(max(rgcca_args$ncomp)))
  }

  func_out$a <- shave(func_out$a, rgcca_args$ncomp)
  func_out$Y <- shave(func_out$Y, rgcca_args$ncomp)

  if (!rgcca_args$superblock) {
    for (j in seq_along(blocks)) {
      rownames(func_out$astar[[j]]) <- colnames(blocks[[j]])
    }
    func_out$astar <- shave(func_out$astar, rgcca_args$ncomp)
  } else {
    rownames(func_out$astar) <- colnames(blocks[[length(blocks)]])
  }

  names(func_out$a) <- names(blocks)
  names(func_out$Y) <- names(blocks)
  if (!rgcca_args$superblock) names(func_out$astar) <- names(blocks)

  is_optimal <- any(rgcca_args[[opt$par]] == "optimal")
  func_out[["optimal"]] <- is_optimal

  if (is_optimal) {
    rgcca_args[[opt$par]] <- func_out$tau
  }

  if (NCOL(rgcca_args[[opt$par]]) > 1) {
    colnames(rgcca_args[[opt$par]]) <- names(blocks)
  }

  if (!is.null(func_out$tau)) {
    func_out$tau <- NULL
  }

  func_out$call <- rgcca_args
  func_out$blocks <- blocks
  func_out$disjunction <- disjunction

  return(func_out)
}
