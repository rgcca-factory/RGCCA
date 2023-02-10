#' Stability selection for SGCCA
#'
#' This function can be used to identify the most stable variables
#' identified as relevant by SGCCA. A Variable Importance in the Projection
#' (VIP) based criterion is used to identify the most stable variables.
#'
#' @inheritParams rgcca_bootstrap
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
#' blocks[[3]] <- Loc
#'
#'
#' fit.sgcca <- rgcca(blocks,
#'   sparsity = c(.071, .2, 1),
#'   ncomp = c(1, 1, 1),
#'   scheme = "centroid",
#'   verbose = TRUE, response = 3
#' )
#'
#' boot.out <- rgcca_bootstrap(fit.sgcca, n_boot = 100, n_cores = 1)
#'
#' fit.stab <- rgcca_stability(fit.sgcca,
#'   keep = sapply(fit.sgcca$a, function(x) mean(x != 0)),
#'   n_cores = 1, n_boot = 10,
#'   verbose = TRUE
#' )
#'
#' boot.out <- rgcca_bootstrap(
#'   fit.stab, n_boot = 500, n_cores = 1, verbose = FALSE
#' )
#' }
#' @export
rgcca_stability <- function(rgcca_res,
                            keep = vapply(
                              rgcca_res$a, function(x) mean(x != 0),
                              FUN.VALUE = 1.0
                            ),
                            n_boot = 100,
                            n_cores = 1,
                            verbose = TRUE,
                            balanced = TRUE,
                            keep_all_variables = FALSE) {
  stopifnot(tolower(rgcca_res$call$method) %in% sparse_methods())
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
      rgcca_bootstrap_k(
        rgcca_res = rgcca_res,
        inds = b, type = "AVE"
      )
    },
    n_cores = n_cores, verbose = verbose
  )

  W <- W[!vapply(W, is.null, logical(1L))]

  res <- format_bootstrap_list(W, rgcca_res)
  J <- length(rgcca_res$blocks)

  if (rgcca_res$call$superblock == TRUE) {
    res <- res[res$block != names(rgcca_res$blocks)[J], ]
    rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[-J]
    rgcca_res$call$blocks <- rgcca_res$call$blocks[-J]
  }

  if (rgcca_res$opt$disjunction) {
     res <- res[res$block != names(rgcca_res$blocks)[rgcca_res$call$response], ]
     rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[-rgcca_res$call$response]
  }

  res_AVE <- res[res$type != "weights", ]
  res <- res[res$type == "weights", ]

  # Compute var2block to later retrieve "block" from "var"
  var2block <- subset(res, res$comp == 1 & res$boot == 1)[, c("var", "block")]
  rownames(var2block) <- var2block$var
  var2block$var <- NULL

  res$scores <- res$value^2 * res_AVE$value
  top <- tapply(
    res$scores, list(var = res$var), mean
  )
  top <- cbind(top = top, block = var2block[names(top), ])

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

  # Keep a percentage of the variables with the top intensities
  keepVar <- lapply(seq_along(rgcca_res$AVE$AVE_X), function(j) {
    x <- top[top[, "block"] == j, "top"]
    order(x, decreasing = TRUE)[seq(round(perc[j] * length(x)))]
  })

  if (rgcca_res$opt$disjunction) {
    keepVar[[rgcca_res$call$response]] <- 1
  }

  rgcca_res$call$blocks <- Map(
    function(x, y) x[, y, drop = FALSE], rgcca_res$call$blocks, keepVar
  )
  rgcca_res$call$tau <-
    rgcca_res$call$sparsity <- rep(1, length(rgcca_res$call$blocks))

  rgcca_res <- rgcca(rgcca_res)

  return(structure(list(
    top = top,
    keepVar = keepVar,
    bootstrap = res,
    rgcca_res = rgcca_res
  ),
  class = "stability"
  ))
}
