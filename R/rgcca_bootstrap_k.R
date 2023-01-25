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
      y <- x[inds, , drop = FALSE]
      rownames(y) <- paste("S", seq_along(inds))
      return(y)
    })
  }
  rgcca_res_boot <- rgcca(rgcca_res)

  # block-weight vector
  missing_var <- unlist(lapply(
    seq_along(rgcca_res_boot$a),
    function(x) {
      setdiff(
        colnames(rgcca_res$blocks[[x]]),
        rownames(rgcca_res_boot$a[[x]])
      )
    }
  ))
  if (length(missing_var) == 0) {
    # block-loadings vector
    A <- check_sign_comp(rgcca_res, rgcca_res_boot$a)

    if (type == "loadings") {
      Y <- lapply(
        seq_along(A),
        function(j) pm(rgcca_res_boot$blocks[[j]], A[[j]])
      )
      L <- lapply(
        seq_along(A),
        function(j) {
          cor(rgcca_res_boot$blocks[[j]], Y[[j]],
              use = "pairwise.complete.obs"
          )
        }
      )
    } else {
      L <- lapply(seq_along(A), function(j) {
        if (j > length(rgcca_res$AVE$AVE_X)) {
          res <- NA
        } else {
          res <-  rgcca_res$AVE$AVE_X[[j]]
        }
        res <- matrix(res, nrow = nrow(A[[j]]), ncol = length(res))
        rownames(res) <- rownames(A[[j]])
        return(res)
      })
    }
    names(L) <- names(rgcca_res$a)
    return(list(W = A, L = L))
  } else {
    return(NULL)
  }
}
