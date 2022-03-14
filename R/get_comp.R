# Get the components of the analysis
#
# @inheritParams plot_ind
# @inheritParams get_ctr
# @param i_blocks A vector of integers giving the index of the blocks
# @param comps A vector of integers giving the index of the components
# @return A matrix containing each selected components and an associated
# response
get_comp <- function(rgcca_res,
                     resp = rep(1, NROW(rgcca_res$Y[[1]])),
                     comps = c(1, 2),
                     i_blocks = rep(length(rgcca_res$Y), 2),
                     predicted = NULL) {
  ### Apply checks
  stopifnot(is(rgcca_res, "rgcca"))
  resp <- as.matrix(check_response(resp))
  if (length(comps) != length(i_blocks)) {
    stop_rgcca(
      "comps and i_blocks must have the same length."
    )
  }

  # Check that i_blocks correspond to existing number of blocks
  lapply(i_blocks, function(i) {
    check_blockx("block index", i, rgcca_res$call$blocks)
    check_ncol(rgcca_res$Y, i)
  })

  # Check that comps and i_blocks correspond to existing components
  mapply(
    function(x, y) check_compx(x, x, rgcca_res$call$ncomp, y), comps, i_blocks
  )

  ### Extract interesting components
  df <- data.frame(
    rgcca_res$Y[[i_blocks[1]]][, comps[1]],
    rgcca_res$Y[[i_blocks[2]]][, comps[2]]
  )

  # If predicted is not null, we overwrite resp
  if (!is.null(predicted)) {
    df2 <- data.frame(
      predicted[[2]][[i_blocks[1]]][, comps[1]],
      predicted[[2]][[i_blocks[2]]][, comps[2]]
    )
    colnames(df2) <- colnames(df)
    df1 <- df[rownames(df2), ]
    df <- rbind(df1, df2)
    df$resp <- as.factor(c(rep("obs", NROW(df2)), rep("pred", NROW(df2))))
    return(df)
  }

  # Otherwise we align resp on Y
  if (is.null(rownames(resp))) {
    if (nrow(resp) != nrow(rgcca_res$Y[[i_blocks[1]]])) {
      stop_rgcca(
        "response must have ", nrow(rgcca_res$Y[[i_blocks[1]]]), " rows."
      )
    }
    df$resp <- resp
  } else {
    df$resp <- apply(
      resp, 2, function(x) x[rownames(rgcca_res$Y[[i_blocks[1]]])]
    )
  }
  return(df)
}
