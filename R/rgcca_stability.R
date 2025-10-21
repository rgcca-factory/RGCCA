#' Identify the most stable variables with SGCCA
#'
#' This function can be used to identify the most stable variables
#' identified as relevant by SGCCA. A Variable Importance in the Projection
#' (VIP) based criterion is used to identify the most stable variables.
#'
#' @inheritParams rgcca_bootstrap
#' @param keep A numeric vector indicating the proportion of variables per
#' block to select.
#' @param verbose A logical value indicating if the progress of the procedure
#' is reported.
#' @return A rgcca_stability object that can be printed and plotted.
#' @return \item{top}{A data.frame giving the indicator (VIP)
#' on which the variables are ranked.}
#' @return \item{n_boot}{The number of bootstrap samples, returned
#' for further use.}
#' @return \item{keepVar}{The indices of the most stable variables.}
#' @return \item{bootstrap}{A data.frame with the block weight vectors
#' computed on each bootstrap sample.}
#' @return \item{rgcca_res}{An RGCCA object fitted on the most stable
#' variables.}
#' @examples
#' \dontrun{
#'  ###########################
#'  # stability and bootstrap #
#'  ###########################
#'
#'  data("ge_cgh_locIGR", package = "gliomaData")
#'  blocks <- ge_cgh_locIGR$multiblocks
#'  Loc <- factor(ge_cgh_locIGR$y)
#'  levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#'  blocks[[3]] <- Loc
#'
#'  fit_sgcca <- rgcca(blocks,
#'     sparsity = c(.071, .2, 1),
#'     ncomp = c(1, 1, 1),
#'     scheme = "centroid",
#'     verbose = TRUE, response = 3
#' )
#'
#'  boot_out <- rgcca_bootstrap(fit_sgcca, n_boot = 100, n_cores = 1)
#'
#'  fit_stab <- rgcca_stability(fit_sgcca,
#'    keep = sapply(fit_sgcca$a, function(x) mean(x != 0)),
#'    n_cores = 1, n_boot = 10,
#'    verbose = TRUE
#'  )
#'
#'  boot_out <- rgcca_bootstrap(
#'    fit_stab, n_boot = 500, n_cores = 1, verbose = TRUE
#'  )
#'
#'  plot(boot_out, block = 1:2, n_mark = 2000, display_order = FALSE)
#' }
#' @export
rgcca_stability <- function(rgcca_res,
                            keep = vapply(
                              rgcca_res$a, function(x) mean(x != 0),
                              FUN.VALUE = 1.0
                            ),
                            n_boot = 100,
                            n_cores = 1,
                            verbose = TRUE) {
  stopifnot(tolower(rgcca_res$call$method) %in% sparse_methods())
  check_integer("n_boot", n_boot)
  check_integer("n_cores", n_cores, min = 0)

  ### Create bootstrap samples
  v_inds <- lapply(seq_len(n_boot), function(i) {
    sample(seq_len(NROW(rgcca_res$call$blocks[[1]])), replace = TRUE)
  })

  ### Run RGCCA on the bootstrap samples
  W <- par_pblapply(v_inds, function(b) {
    rgcca_bootstrap_k(
      rgcca_res = rgcca_res,
      inds = b, type = "AVE"
    )
  }, n_cores = n_cores, verbose = verbose)

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
    n_boot = n_boot,
    keepVar = keepVar,
    bootstrap = res,
    rgcca_res = rgcca_res
  ),
  class = "rgcca_stability"
  ))
}
