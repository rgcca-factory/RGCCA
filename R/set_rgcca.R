#' Set and run RGCCA on a subset of individuals (used in bootstrap and
#' cross validation procedures)
#'
#' set_rgcca() runs a new RGCCA model with the same parameters as a fitted RGCCA
#' model (or other ones if specified) on a subset of individuals.
#' @param inds vector of individuals removed from the blocks
#' @param keep_inds logical, if TRUE, inds are kept instead of being removed
#' @noRd
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
  ### Copy parameters from rgcca_res$call if not overwritten, parameters
  #   blocks, response and keep_inds are treated differently
  all_args <- names(environment())
  used_args <- c(names(match.call()), "blocks", "response", "keep_inds")
  for (n in setdiff(all_args, used_args)) {
    assign(n, rgcca_res$call[[n]])
  }

  if (is.null(blocks)) {
    blocks <- rgcca_res$call$raw

    if (!is.null(rgcca_res$call$response)) {
      response <- rgcca_res$call$response
    }
  }

  method <- rgcca_res$call$method

  ### Subset blocks to remove inds
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

  ### Call rgcca
  return(rgcca(
    boot_blocks,
    connection = connection,
    tau = tau,
    sparsity = sparsity,
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
  ))
}
