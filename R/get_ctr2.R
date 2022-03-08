# Get the indexes of the analysis
#
# @inheritParams plot_var_2D
# @inheritParams get_ctr
# @return A matrix containg the indexes (correlation of the blocks with a
# component or their weights) for each selected component and an associated
# response

get_ctr2 <- function(rgcca_res,
                     compx = 1,
                     compy = 2,
                     compz = NULL,
                     i_block = length(rgcca_res$call$blocks),
                     type = "loadings",
                     n_mark = 100,
                     collapse = FALSE,
                     remove_var = TRUE,
                     resp = NULL) {
  stopifnot(is(rgcca_res, "rgcca"))
  check_blockx("i_block", i_block, rgcca_res$call$blocks)
  check_ncol(rgcca_res$a, i_block)
  for (i in c("compx", "compy", "compz")) {
    if (!is.null(get(i))) {
      check_compx(i, get(i), rgcca_res$call$ncomp, i_block)
    }
  }
  check_integer("n_mark", n_mark)
  for (i in c("collapse", "remove_var")) {
    check_boolean(i, get(i))
  }

  selected_var <- NULL
  blocks <- rgcca_res$call$blocks

  if (collapse) {
    if (rgcca_res$call$superblock) {
      blocks <- blocks[-length(blocks), drop = FALSE]
      if (i_block > length(blocks)) {
        i_block <- length(blocks)
      }
    }
    blocks_all <- blocks
    blocks <- rep(list(Reduce(cbind, blocks)), length(blocks))
    names(blocks) <- names(blocks_all)
  }

  df <- get_ctr(rgcca_res, compx, compy, compz, i_block, type, collapse)

  if (tolower(rgcca_res$call$method) %in% c("spls", "spca", "sgcca")) {
    if (collapse) {
      J <- seq(length(rgcca_res$a))
    } else {
      J <- i_block
    }

    selected_var <- unlist(
      lapply(
        J,
        function(x) {
          apply(
            sapply(
              c(compx, compy, compz[compz >= rgcca_res$call$ncomp[x]]),
              function(y) rgcca_res$a[[x]][, y] != 0
            ),
            1,
            function(z) Reduce("|", z)
          )
        }
      )
    )
    df <- df[names(which(selected_var)), ]
  }

  if (n_mark > NROW(df)) {
    n_mark <- NROW(df)
  }

  # TODO: function in other place
  if (remove_var) {
    selected_var <- unique(as.vector(
      sapply(seq(length(c(compx, compy, compz))), function(x) {
        row.names(data.frame(
          df[order(abs(df[, x]), decreasing = TRUE), ]
        )[seq(n_mark), ])
      })
    ))
    df <- df[selected_var, ]
  } else {
    selected_var <- row.names(df)
  }

  # group by blocks
  if (is.null(resp)) {
    if (
      (rgcca_res$call$superblock && i_block == length(rgcca_res$a)) ||
        collapse
    ) {
      if (collapse) {
        resp <- get_bloc_var(lapply(blocks_all, t), TRUE)
      } else {
        resp <- get_bloc_var(rgcca_res$a)

        resp <- resp[
          unlist(
            lapply(
              seq(length(selected_var)),
              function(x) {
                which(
                  colnames(blocks[[length(blocks)]]) == selected_var[x]
                )
              }
            )
          )
        ]
      }
      resp <- resp[row.names(df)]
    } else {
      resp <- rep(1, NROW(df))
    }
  }

  data.frame(df, resp)
}
