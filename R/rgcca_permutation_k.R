#' An internal function used by rgcca_permutation() to perform multiple s/rgcca
#' on permuted blocks
#'
#' If superblock is TRUE, blocks are permuted and then a superblock is created
#' by concatenating the permuted blocks.
#' @inheritParams rgcca_permutation
#' @noRd
rgcca_permutation_k <- function(rgcca_args, perm, par_type, par_value) {
  if (perm) {
    rgcca_args$blocks <- lapply(rgcca_args$blocks, function(x) {
      block <- as.matrix(x)[sample(seq_len(NROW(x))), , drop = FALSE]
      rownames(block) <- rownames(x)
      return(block)
    })
  }

  rgcca_args[[par_type]] <- par_value

  res <- do.call(rgcca, rgcca_args)

  if (max(res$call$ncomp) > 1) {
    criterion <- vapply(res$crit, function(x) {
      x[length(x)]
    }, FUN.VALUE = numeric(1))
    crit_permut <- sum(criterion)
  } else {
    crit_permut <- res$crit[length(res$crit)]
  }

  if (perm) {
    par_value <- c()
  } else {
    par_value <- res$call[[par_type]]
  }

  return(list(crit = crit_permut, par_value = par_value))
}
