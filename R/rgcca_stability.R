#' Stability selection for SGCCA
#'
#' This function can be used to identify the most stable variables
#' identified as relevant by SGCCA. A Variable Importance in the Projection
#' (VIP) based criterion is used to identify the most stable variables.
#'
#' @inheritParams bootstrap
#' @param rgcca_res A fitted RGCCA object (see \code{\link[RGCCA]{rgcca}})
#' @param keep numeric vector indicating the proportion of top variables per
#' block.
#' @param n_boot Number of bootstrap samples (Default: 100).
#' @param n_cores Number of cores for parallelization.
#' @param verbose Logical value indicating if the progress of the procedure
#' is reported.
#' @return \item{top}{indicator on which variables are ranked.}
#' @return \item{keepVar}{indices of the top variables.}
#' @return \item{bootstrap}{block-weight vectors for ech bootstrap sample.}
#' @return \item{rgcca_res}{an RGCCA object fitted on the most stable
#' variables.}
#' @examples
#' \dontrun{
#' ###########################
#' # stability and bootstrap #
#' ###########################
#'
#' data("ge_cgh_locIGR", package = "gliomaData")
#' blocks <- ge_cgh_locIGR$multiblocks
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' connection <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' blocks[[3]] <- blocks[[3]][, -3]
#'
#' fit.sgcca <- rgcca(blocks,
#'   connection = connection,
#'   sparsity = c(.071, .2, 1),
#'   ncomp = c(1, 1, 1),
#'   scheme = "centroid",
#'   verbose = TRUE
#' )
#'
#' boot.out <- bootstrap(fit.sgcca, n_boot = 100, n_cores = 1, verbose = TRUE)
#'
#' fit.stab <- rgcca_stability(fit.sgcca,
#'   keep = sapply(fit.sgcca$a, function(x) mean(x != 0)),
#'   n_cores = 15,
#'   verbose = TRUE
#' )
#'
#' boot.out <- bootstrap(fit.stab, n_boot = 500, n_cores = 1, verbose = FALSE)
#' }
#' @export
rgcca_stability <- function(rgcca_res,
                            keep = vapply(
                              rgcca_res$a, function(x) mean(x != 0),
                              FUN.VALUE = 1.0
                            ),
                            n_boot = 100,
                            n_cores = 1,
                            verbose = FALSE,
                            balanced = TRUE,
                            keep_all_variables = FALSE) {
  stopifnot(tolower(rgcca_res$call$method) %in% c("sgcca", "spls", "spca"))
  check_integer("n_boot", n_boot)
  check_integer("n_cores", n_cores, min = 0)

  boot_sampling <- generate_resampling(
    rgcca_res = rgcca_res,
    n_boot = n_boot,
    balanced = balanced,
    verbose = verbose,
    keep_all_variables = keep_all_variables
  )

  sd_null <- boot_sampling$sd_null

  if (!is.null(sd_null)) {
    rgcca_res$call$blocks <- remove_null_sd(
      list_m = rgcca_res$call$blocks,
      column_sd_null = sd_null
    )$list_m
    rgcca_res <- rgcca(rgcca_res)
  }

  W <- par_pblapply(
    boot_sampling$full_idx, function(b) {
      bootstrap_k(
        rgcca_res = rgcca_res,
        inds = b
      )
    },
    n_cores = n_cores, verbose = verbose
  )

  list_res <- format_bootstrap_list(W, rgcca_res, n_boot, 1)

  J <- length(list_res[[1]])

  if (rgcca_res$call$superblock == TRUE) {
    list_res <- lapply(list_res, function(x) x[-length(x)])
    rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[-J]
    rgcca_res$call$blocks <- rgcca_res$call$blocks[-J]
  }

  if (isTRUE(rgcca_res$call$disjunction)) {
     list_res <- lapply(list_res, function(x) x[-rgcca_res$call$response])
     rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[-rgcca_res$call$response]
  }

  mylist <- lapply(
    seq_along(list_res),
    function(i) {
      lapply(
        list_res[[i]],
        function(x) {
          apply(x, 1, function(y) sum(abs(y), na.rm = TRUE))
        }
      )
    }
  )

  intensity <- Map(
    "*",
    do.call(Map, c(rbind, mylist)),
    rgcca_res$AVE$AVE_X
  )

  top <- lapply(intensity, function(x) colMeans(x, na.rm = TRUE))
  perc <- elongate_arg(keep, top)

  if (is.null(dim(rgcca_res$call$sparsity))) {
    if (rgcca_res$call$superblock == TRUE) {
      rgcca_res$call$sparsity <- rgcca_res$call$sparsity[-J]
    }
    perc[which(rgcca_res$call$sparsity == 1)] <- 1
  } else {
    if (rgcca_res$call$superblock == TRUE) {
      rgcca_res$call$sparsity <- rgcca_res$call$sparsity[, -J]
    }
    perc[which(rgcca_res$call$sparsity[1, ] == 1)] <- 1
  }

  keepVar <- lapply(
    seq_along(top),
    function(x) {
      order(top[[x]],
        decreasing = TRUE
      )[1:round(perc[x] * length(top[[x]]))]
    }
  )

  if (isTRUE(rgcca_res$call$disjunction)) {
    keepVar[[rgcca_res$call$response]] <- 1
  }

  rgcca_res$call$blocks <- Map(
    function(x, y) x[, y], rgcca_res$call$blocks, keepVar
  )
  rgcca_res$call$tau <-
    rgcca_res$call$sparsity <- rep(1, length(rgcca_res$call$blocks))

  rgcca_res <- rgcca(rgcca_res)

  return(structure(list(
    top = top,
    keepVar = keepVar,
    bootstrap = list_res,
    rgcca_res = rgcca_res
  ),
  class = "stability"
  ))
}
