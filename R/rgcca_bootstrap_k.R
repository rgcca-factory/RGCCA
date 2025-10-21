#' Internal function for computing bootstrap of RGCCA.
#' @inheritParams rgcca_bootstrap
#' @param inds A vector of integers defining the index of the observations
#' taken into account for this bootstrap sample.
#' @param type A character indicating the type of the second object to return.
#' @return \item{W}{A list of RGCCA bootstrap weights. Returned only if there
#' are no missing variables, otherwise the name of the missing variables are
#' returned.}
#' @return \item{L}{If type == "loadings", a list of RGCCA bootstrap loadings.
#' Returned only if there are no missing variables, otherwise the name of the
#' missing variables are returned.
#' If type == "AVE", the AVE of the fitted RGCCA model.}
#' @title Compute bootstrap (internal).
#' @noRd
rgcca_bootstrap_k <- function(rgcca_res, inds = NULL, type = "loadings") {
  if (length(inds) > 0) {
    rgcca_res$call$blocks <- lapply(rgcca_res$call$blocks, function(x) {
      y <- subset_block_rows(x, inds, drop = FALSE)
      rownames(y) <- paste("S", seq_along(inds))
      return(y)
    })
  }
  rgcca_res_boot <- rgcca(rgcca_res)

  # block-loadings vector
  A <- check_sign_comp(rgcca_res, rgcca_res_boot$a)

  if (type == "loadings") {
    Y <- lapply(
      seq_along(A),
      function(j) pm(to_mat(rgcca_res_boot$blocks[[j]]), A[[j]])
    )
    L <- lapply(
      seq_along(A),
      function(j) {
        cor2(to_mat(rgcca_res_boot$blocks[[j]]), Y[[j]])
      }
    )
  } else {
    L <- lapply(names(A), function(n) {
      if (!(n %in% names(rgcca_res$AVE$AVE_X))) {
        res <- rep(-1, NCOL(rgcca_res$a[[n]]))
      } else {
        res <-  rgcca_res_boot$AVE$AVE_X[[n]]
      }
      res <- matrix(
        res, nrow = nrow(A[[n]]),
        ncol = length(res), byrow = TRUE
      )
      rownames(res) <- rownames(A[[n]])
      return(res)
    })
  }
  names(L) <- names(rgcca_res$a)
  return(list(W = A, L = L))
}
