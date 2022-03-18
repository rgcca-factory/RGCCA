#' Stability selection for SGCCA
#'
#' This function can be used to identify the most stable variables
#' identified as relevant by SGCCA. The Variable Importance in the Projection
#' (VIP) criterion is used to identify the most stable variables.
#'
#' @param rgcca_res A fitted RGCCA object (see \code{\link[RGCCA]{rgcca}})
#' @param keep numeric vector indicating the proportion of top variables per
#' block.
#' @param n_boot Number of bootstrap samples (Default: 100).
#' @param n_cores Number of cores for parallelization.
#' @param verbose Logical value indicating if the progress of the procedure
#' is reported.
#' @inheritParams generate_resampling
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
#' require(gliomaData)
#' data(ge_cgh_locIGR)
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
                            keep = sapply(
                              rgcca_res$a, function(x) mean(x != 0)
                            ),
                            n_boot = 100,
                            n_cores = parallel::detectCores() - 1,
                            verbose = FALSE,
                            balanced = TRUE,
                            keep_all_variables = FALSE) {
  local_env <- new.env()

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
    rgcca_res$call$raw <- remove_null_sd(
      list_m = rgcca_res$call$raw,
      column_sd_null = sd_null
    )$list_m
    rgcca_res <- set_rgcca(rgcca_res)
  }

  ndefl_max <- max(rgcca_res$call$ncomp)
  list_res <- list()
  for (i in 1:ndefl_max) {
    list_res[[i]] <- list()
    for (block in names(rgcca_res$call$blocks)) {
      list_res[[i]][[block]] <-
        matrix(NA, dim(rgcca_res$call$blocks[[block]])[2], n_boot)
      rownames(list_res[[i]][[block]]) <-
        colnames(rgcca_res$call$blocks[[block]])
    }
  }

  if (n_cores == 0) n_cores <- 1

  if (!verbose) {
    pbapply::pboptions(type = "none")
  } else {
    pbapply::pboptions(type = "timer")
  }

  if (Sys.info()["sysname"] == "Windows") {
    if (n_cores > 1) {
      assign("rgcca_res", rgcca_res, envir = local_env)
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, "rgcca_res", envir = local_env)
      W <- pbapply::pblapply(boot_sampling$full_idx,
        function(b) {
          bootstrap_k(
            rgcca_res = rgcca_res,
            inds = b
          )
        },
        cl = cl
      )
      parallel::stopCluster(cl)
    } else {
      W <- pbapply::pblapply(
        boot_sampling$full_idx,
        function(b) {
          bootstrap_k(rgcca_res = rgcca_res, inds = b)
        }
      )
    }
  } else {
    W <- pbapply::pblapply(boot_sampling$full_idx,
      function(b) {
        bootstrap_k(rgcca_res = rgcca_res, inds = b)
      },
      cl = n_cores
    )
  }

  W <- lapply(W, `[[`, 1)

  for (k in seq(n_boot)) {
    for (i in seq(ndefl_max)) {
      for (j in seq_along(rgcca_res$call$blocks)) {
        block <- names(rgcca_res$call$blocks)[j]
        if (!is.null(names(W[[k]]))) {
          if (i <= NCOL(W[[k]][[block]])) {
            list_res[[i]][[block]][, k] <- W[[k]][[block]][, i]
          }
        } else {
          list_res[[i]][[block]][, k] <-
            rep(NA, length(list_res[[i]][[block]][, k]))
          if (is.character(W[[k]])) {
            warning(paste0(
              "This bootstrap sample was removed due to zero variance ",
              "variable(s): ",
              paste(W[[k]], collapse = " - ")
            ))
          }
        }
      }
    }
  }

  J <- length(list_res[[1]])

  if (rgcca_res$call$superblock == TRUE) {
    list_res <- lapply(list_res, function(x) x[-length(x)])
    rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[-J]
    rgcca_res$call$raw <- rgcca_res$call$raw[-J]
  }

  mylist <- lapply(
    seq_along(list_res),
    function(i) {
      sapply(
        list_res[[i]],
        function(x) {
          apply(x, 1, function(y) sum(abs(y), na.rm = TRUE))
        }
      )
    }
  )

  intensity <- mapply(
    "*",
    do.call(mapply, c(rbind, mylist)),
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

  new_block <- mapply(function(x, y) x[, y], rgcca_res$call$raw, keepVar,
    SIMPLIFY = FALSE
  )

  rgcca_res <- rgcca(new_block,
    connection = rgcca_res$call$connection,
    superblock = rgcca_res$call$superblock,
    ncomp = rgcca_res$call$ncomp,
    bias = rgcca_res$call$bias,
    tau = 1,
    scale = rgcca_res$call$scale,
    verbose = FALSE,
    scale_block = rgcca_res$call$scale_block
  )

  return(structure(list(
    top = top,
    keepVar = keepVar,
    bootstrap = list_res,
    rgcca_res = rgcca_res
  ),
  class = "stability"
  ))
}
