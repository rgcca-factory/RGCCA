# inds : individuals removed from the blocks (for crossvalidation)
# affects same parameters as rgcca_res (or other ones if specified) on a subset
# of individuals determined by inds (individuals to remove)
set_rgcca <- function(rgcca_res,
                      blocks = NULL,
                      connection = NULL,
                      tau = NULL,
                      sparsity = NULL,
                      ncomp = NULL,
                      scheme = NULL,
                      init = NULL,
                      bias = TRUE,
                      tol = 1e-03,
                      scale = NULL,
                      scale_block = NULL,
                      superblock = NULL,
                      response = NULL,
                      NA_method = NULL,
                      inds = NULL,
                      keep_inds = FALSE) {
  if (is.null(connection)) {
    connection <- rgcca_res$call$connection
  }
  if (is.null(scale)) {
    scale <- rgcca_res$call$scale
  }
  if (is.null(scale_block)) {
    scale_block <- rgcca_res$call$scale_block
  }
  if (is.null(superblock)) {
    superblock <- rgcca_res$call$superblock
  }
  if (is.null(NA_method)) {
    NA_method <- rgcca_res$call$NA_method
  }
  if (is.null(scheme)) {
    scheme <- rgcca_res$call$scheme
  }
  if (is.null(bias)) {
    bias <- rgcca_res$call$bias
  }
  if (is.null(init)) {
    init <- rgcca_res$call$init
  }
  if (is.null(ncomp)) {
    ncomp <- rgcca_res$call$ncomp
  }
  if (is.null(sparsity)) {
    sparsity <- rgcca_res$call$sparsity
  }
  if (is.null(tau)) {
    tau <- rgcca_res$call$tau
  }
  if (is.null(blocks)) {
    blocks <- rgcca_res$call$raw

    if (!is.null(rgcca_res$call$response)) {
      response <- rgcca_res$call$response
    }
  }

  method <- rgcca_res$call$method
  if (tolower(method) %in% c("sgcca", "spca", "spls")) {
    par <- "sparsity"
    penalty <- sparsity
  } else {
    par <- "tau"
    penalty <- tau
  }

  if (length(inds) == 0) {
    boot_blocks <- blocks
  } else {
    if (keep_inds) {
      boot_blocks <- lapply(blocks, function(x) {
        y <- x[inds, , drop = FALSE]
        rownames(y) <- paste("S", seq_along(inds))
        return(y)
      })
    } else {
      boot_blocks <- lapply(blocks, function(x) x[-inds, , drop = FALSE])
      if ("character" %in% class(boot_blocks[[response]])) {
        if (length(unique(boot_blocks[[response]])) == 1) {
          warning("One block has no variablity and rgcca fails to fit.")
          return(NULL)
        }
      }
    }
  }

  func <- quote(
    rgcca(
      boot_blocks,
      connection = connection,
      superblock = superblock,
      response = response,
      ncomp = ncomp,
      scheme = scheme,
      scale = scale,
      scale_block = scale_block,
      method = method,
      verbose = FALSE,
      init = init,
      bias = bias,
      NA_method = NA_method,
      tol = tol
    )
  )

  func[[par]] <- penalty

  res <- eval(as.call(func))

  return(res)
}
