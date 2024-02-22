#' Format output of rgcca function
#'
#' @noRd
format_output <- function(func_out, rgcca_args, opt, blocks) {
  ### Compute AVE
  blocks_AVE <- seq_along(blocks)
  names_AVE <- names(blocks)
  ncomp_AVE <- rgcca_args$ncomp
  if (opt$disjunction) {
    blocks_AVE <- blocks_AVE[-rgcca_args$response]
    names_AVE <- names_AVE[-rgcca_args$response]
    ncomp_AVE <- ncomp_AVE[-rgcca_args$response]
  }

  # Compute the proportion of variance explained between blocks
  # If a component is null (all elements are 0), the correlation
  # with the other blocks is set to 0
  AVE_inner <- vapply(seq(max(rgcca_args$ncomp)), function(n) {
    cor_matrix <- cor2(
      do.call(cbind, lapply(func_out$Y, function(y) y[, n]))
    )
    sum(rgcca_args$connection * cor_matrix^2) / (sum(rgcca_args$connection))
  }, FUN.VALUE = double(1L))

  AVE <- lapply(blocks_AVE, function(j) {
    ave(blocks[[j]], func_out$Y[[j]])
  })
  AVE_X <- lapply(AVE, "[[", "AVE_X")
  AVE_X_cor <- lapply(AVE, "[[", "AVE_X_cor")
  var_tot <- vapply(AVE, "[[", "var_tot", FUN.VALUE = 1.)

  outer <- matrix(unlist(AVE_X_cor), nrow = max(ncomp_AVE))
  AVE_outer <- as.vector((outer %*% var_tot) / sum(var_tot))
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

  if (rgcca_args$superblock && rgcca_args$comp_orth) {
    rownames(func_out$astar) <- colnames(blocks[[length(blocks)]])
  } else {
    for (j in seq_along(blocks)) {
      rownames(func_out$astar[[j]]) <- colnames(blocks[[j]])
    }
    func_out$astar <- shave(func_out$astar, rgcca_args$ncomp)
    names(func_out$astar) <- names(blocks)
  }

  names(func_out$a) <- names(blocks)
  names(func_out$Y) <- names(blocks)

  is_optimal <- any(rgcca_args[[opt$param]] == "optimal")
  func_out[["optimal"]] <- is_optimal

  if (is_optimal) {
    rgcca_args[[opt$param]] <- func_out$tau
  }

  if (NCOL(rgcca_args[[opt$param]]) > 1) {
    colnames(rgcca_args[[opt$param]]) <- names(blocks)
  }

  if (!is.null(func_out$tau)) {
    func_out$tau <- NULL
  }

  func_out$opt <- opt
  func_out$call <- rgcca_args
  func_out$blocks <- blocks

  return(func_out)
}
