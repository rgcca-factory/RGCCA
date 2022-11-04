#' An internal function used by rgcca_permutation() to perform multiple s/rgcca
#' on permuted blocks
#'
#' If superblock is TRUE, blocks are permuted and then a superblock is created
#' by concatenating the permuted blocks.
#' @inheritParams rgcca_permutation
#' @noRd
rgcca_permutation_k <- function(blocks,
                                method = "rgcca",
                                scale = TRUE,
                                scale_block = TRUE,
                                connection = 1 - diag(length(blocks)),
                                scheme = "factorial",
                                ncomp = rep(1, length(blocks)),
                                tau = rep(1, length(blocks)),
                                sparsity = rep(1, length(blocks)),
                                init = "svd",
                                bias = TRUE,
                                tol = 1e-08,
                                response = NULL,
                                superblock = FALSE,
                                NA_method = "nipals",
                                quiet = TRUE,
                                perm = TRUE,
                                par_type = "tau",
                                par_value = rep(1, length(blocks))) {
  blocks_to_use <- blocks
  if (perm) {
    blocks_to_use <- lapply(
      blocks,
      function(x) {
        blocks_to_use_k <-
          as.matrix(x)[sample(seq_len(NROW(x))), , drop = FALSE]
        rownames(blocks_to_use_k) <- rownames(x)
        return(blocks_to_use_k)
      }
    )
  }

  rgcca_args <- list(
    blocks = blocks_to_use,
    method = method,
    scale = scale,
    scale_block = scale_block,
    connection = connection,
    scheme = scheme,
    init = init,
    bias = bias,
    tol = tol,
    response = response,
    superblock = superblock,
    NA_method = NA_method,
    quiet = quiet,
    ncomp = ncomp,
    tau = tau,
    sparsity = sparsity
  )
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

  return(crit_permut)
}
