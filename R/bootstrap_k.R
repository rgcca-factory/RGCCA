#' Internal function for computing bootstrap of RGCCA.
#' @inheritParams bootstrap
#' @param inds A vector of integers defining the index of the observations
#' taken into account for this bootstrap sample.
#' @return \item{W}{A list of RGCCA bootstrap weights. Returned only if there
#' are no missing variables, otherwise the name of the missing variables are
#' returned.}
#' @return \item{L}{A list of RGCCA bootstrap loadings.  Returned only if there
#' are no missing variables, otherwise the name of the missing variables are
#' returned.}
#' @title Compute bootstrap (internal).
#' @noRd
bootstrap_k <- function(rgcca_res, inds = NULL) {
  rgcca_res_boot <- set_rgcca(rgcca_res,
    inds      = inds,
    keep_inds = TRUE,
    NA_method = "nipals"
  )
  # block-weight vector
  missing_var <- unlist(lapply(
    seq_along(rgcca_res_boot$a),
    function(x) {
      setdiff(
        colnames(rgcca_res$call$blocks[[x]]),
        rownames(rgcca_res_boot$a[[x]])
      )
    }
  ))
  if (length(missing_var) == 0) {
    # block-loadings vector
    A <- check_sign_comp(rgcca_res, rgcca_res_boot$a)

    Y <- lapply(
      seq_along(A),
      function(j) pm(rgcca_res_boot$call$blocks[[j]], A[[j]])
    )
    L <- lapply(
      seq_along(A),
      function(j) {
        cor(rgcca_res_boot$call$blocks[[j]], Y[[j]],
          use = "pairwise.complete.obs"
        )
      }
    )

    names(L) <- names(rgcca_res$a)
    return(list(W = A, L = L))
  } else {
    return(list(W = missing_var, L = missing_var))
  }
}
