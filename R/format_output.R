#' Format output of rgcca function
#'
#' @noRd
format_output <- function(func_out, opt, raw, func_call = NULL) {
  ### Compute AVE
  blocks_AVE <- seq_along(opt$blocks)
  names_AVE <- names(opt$blocks)
  ncomp_AVE <- opt$ncomp
  pjs <- vapply(opt$blocks, NCOL, FUN.VALUE = 1L)
  if (!is.null(func_call$disjunction)) {
    blocks_AVE <- blocks_AVE[-opt$response]
    names_AVE <- names_AVE[-opt$response]
    ncomp_AVE <- ncomp_AVE[-opt$response]
    pjs <- pjs[-opt$response]
  }

  AVE_X <- lapply(blocks_AVE, function(j) apply(
    func_out$Y[[j]], 2, rsq, opt$blocks[[j]]
  ))
  AVE_X_cum <- lapply(blocks_AVE, function(j) {
    vapply(
      seq(NCOL(func_out$Y[[j]])),
      function(p) rsq(func_out$Y[[j]][, seq(p)], opt$blocks[[j]]),
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
    AVE_outer = AVE_outer, AVE_inner = func_out$AVE_inner
  )
  func_out$AVE_inner <- NULL

  names(func_out$AVE$AVE_X) <- names_AVE
  names(func_out$AVE$AVE_X_cor) <- names_AVE

  ### Set names and shave
  for (j in seq_along(opt$blocks)) {
    rownames(func_out$a[[j]]) <- colnames(opt$blocks[[j]])
    rownames(func_out$Y[[j]]) <- rownames(opt$blocks[[j]])
    colnames(func_out$Y[[j]]) <- paste0("comp", seq_len(max(opt$ncomp)))
  }

  func_out$a <- shave(func_out$a, opt$ncomp)
  func_out$Y <- shave(func_out$Y, opt$ncomp)

  if (!opt$superblock) {
    for (j in seq(length(opt$blocks))) {
      rownames(func_out$astar[[j]]) <- colnames(opt$blocks[[j]])
    }
    func_out$astar <- shave(func_out$astar, opt$ncomp)
  } else {
    rownames(func_out$astar) <- colnames(opt$blocks[[length(opt$blocks)]])
  }

  names(func_out$a) <- names(opt$blocks)
  names(func_out$Y) <- names(opt$blocks)
  if (!opt$superblock) names(func_out$astar) <- names(opt$blocks)

  func_out$call <- list(
    blocks = opt$blocks,
    connection = opt$connection,
    superblock = opt$superblock,
    ncomp = opt$ncomp,
    scheme = opt$scheme,
    response = opt$response,
    raw = raw,
    method = opt$method
  )

  is_optimal <- any(opt$penalty == "optimal")
  func_out$call[["optimal"]] <- is_optimal

  if (is_optimal) {
    func_out$call[[opt$par]] <- func_out$tau
  } else {
    func_out$call[[opt$par]] <- opt$penalty
  }

  if (NCOL(func_out$call[[opt$par]]) > 1) {
    colnames(func_out$call[[opt$par]]) <- names(opt$blocks)
  }

  if (!is.null(func_out$tau)) {
    func_out$tau <- NULL
  }

  func_out$call <- c(func_out$call, func_call)

  return(func_out)
}
