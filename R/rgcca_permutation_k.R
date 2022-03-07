# An internal function used by rgcca_permutation() to perform multiple s/rgcca
# on permuted blocks
#
# data("Russett")
# blocks = list(agriculture = Russett[, seq(3)],
#               industry = Russett[, 4:5],
#               politic = Russett[, 6:11])
#
# rgcca_permutation_k(blocks)
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
                                rgcca_res = NULL,
                                par_type = "tau",
                                par_value = rep(1, length(blocks))) {
  if (!is.null(rgcca_res)) {
    stopifnot(is(rgcca_res, "rgcca"))
    method <- rgcca_res$call$method
    scale_block <- rgcca_res$call$scale_block
    scale <- rgcca_res$call$scale
    scheme <- rgcca_res$call$scheme
    response <- rgcca_res$call$response
    tol <- rgcca_res$call$tol
    NA_method <- rgcca_res$call$NA_method
    init <- rgcca_res$call$init
    bias <- rgcca_res$call$bias
    blocks <- rgcca_res$call$raw
    superblock <- rgcca_res$call$superblock
    connection <- rgcca_res$call$connection
    tau <- rgcca_res$call$tau
    ncomp <- rgcca_res$call$ncomp
    sparsity <- rgcca_res$call$sparsity
  }

  if (method %in% c("sgcca", "spls")) {
    par_type <- "sparsity"
  } else {
    par_type <- "tau"
  }

  if (perm) {
    blocks_to_use <- blocks
    blocks_to_use <- lapply(
      seq(length(blocks)),
      function(k) {
        blocks_to_use_k <-
          as.matrix(blocks[[k]][sample(seq(NROW(blocks[[k]]))), ])
        rownames(blocks_to_use_k) <- rownames(blocks[[k]])
        return(blocks_to_use_k)
      }
    )
    names(blocks_to_use) <- names(blocks)
  } else {
    blocks_to_use <- blocks
  }

  func <- quote(
    rgcca(
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
      quiet = quiet
    )
  )

  func$ncomp <- ncomp
  func$tau <- tau
  func$sparsity <- sparsity
  func[[par_type]] <- par_value

  res <- eval(as.call(func))

  if (max(ncomp) > 1) {
    criterion <- sapply(res$crit, function(x) {
      x[length(x)]
    })
    crit_permut <- sum(criterion)
  } else {
    crit_permut <- res$crit[length(res$crit)]
  }

  return(crit_permut)
}
