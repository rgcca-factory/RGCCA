# Add removed variables in one submodel (bootstrap, cross-validation)
# to a list of weights and affect them to 0
add_variables_attr <- function(rgcca_res, w, type = "scale") {
  blocks_all <- list()
  for (x in seq_along(rgcca_res$call$blocks)) {
    blocks_all[[x]] <- rgcca_res$call$blocks[[x]][intersect(rownames(rgcca_res$call$blocks[[1]]), rownames(w[[1]])), ]
  }

  blocks_all2 <- lapply(
    1:length(blocks_all),
    function(i) {
      if (is.null(dim(blocks_all[[i]]))) {
        blocks_all[[i]] <- matrix(blocks_all[[i]], ncol = 1)
        colnames(blocks_all[[i]]) <- names(blocks_all)[i]
      }
      if (is.null(colnames(blocks_all[[i]]))) {
        colnames(blocks_all[[i]]) <- names(blocks_all)[i]
      }
      return(blocks_all[[i]])
    }
  )
  names(blocks_all2) <- names(blocks_all)

  blocks_all <- blocks_all2

  missing_var <- lapply(
    seq_along(w),
    function(x) {
      if (!is.null(dim(blocks_all[[x]]))) {
        setdiff(
          colnames(blocks_all[[x]]),
          colnames(w[[x]])
        )
      }
    }
  )

  missing_tab <- lapply(
    seq_along(w),
    function(x) {
      # setNames(

      for (i in seq_along(missing_var[[x]]))
      {
        if (type == "scale") {
          0
        } else
        if (!is.null(dim(blocks_all))) {
          mean(blocks_all[[x]][, i])
        }
      }
    }
  )

  for (i in seq_along(w)) {
    if (NROW(missing_tab[[i]]) != 0) {
      names(missing_tab[[i]]) <- missing_var[[i]]
      w[[i]] <- c(w[[i]], missing_tab[[i]])
    }
  }

  w <- lapply(seq_along(w), function(x) {
    if (!is.null(colnames(blocks_all[[x]]))) {
      w[[x]][colnames(blocks_all[[x]]), drop = FALSE]
    } else {
      w[[x]]
    }
  })


  # names(w) <- names(rgcca_res$call$blocks)
  return(w)
}
