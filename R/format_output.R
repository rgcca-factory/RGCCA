format_output <- function(func_out, opt, raw, func_call = NULL) {
  for (j in seq(length(opt$blocks))) {
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

  names(func_out$AVE$AVE_X) <- names(opt$blocks)

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

  if (!is.null(func_out$tau)) {
    func_out$tau <- NULL
  }

  func_out$call <- c(func_out$call, func_call)

  return(func_out)
}
