#' Bootstrap confidence intervals and p-values
#'
#' Boostrap confidence intervals and p-values for evaluating the significancy/
#' stability of the block-weight vectors produce by S/RGCCA.
#' @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})
#' @param n_boot Number of bootstrap samples. Default: 100.
#' @param n_cores Number of cores for parallelization.
#' @inheritParams generate_resampling
#' @param verbose Logical value indicating if the progress of the
#' bootstrap procedure is reported.
#' @return A list containing two objects: 'bootstrap' and 'rgcca'.
#' 'bootstrap' is a list containing for each block, a matrix
#' with the variables of the block in row and the block weight vector
#' calculated accross bootstrap sample in column. 'rgcca' is the fitted rgcca
#' object obtained from the original data. (see  \code{\link[RGCCA]{rgcca}})
#' @examples
#' # Bootstrap confidence intervals and p-values for RGCCA
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#'
#' fit.rgcca <- rgcca(blocks, ncomp = c(2, 1, 2))
#' boot.out <- bootstrap(fit.rgcca, n_boot = 20, n_cores = 2)
#'
#' plot(boot.out, type = "weight", block = 3, comp = 1)
#'
#' print(boot.out, comp = 2)
#' get_bootstrap(boot.out, block = 1, comp = 1)
#'
#' fit.rgcca <- rgcca(blocks, method = "mcoa")
#' boot.out <- bootstrap(fit.rgcca, n_boot = 50, n_cores = 2)
#'
#' plot(boot.out, type = "weight", block = 4)
#' \dontrun{
#' # Stability of the selected variables for SGCCA
#' # Not run:
#' # Download the dataset's package at http://biodev.cea.fr/sgcca/.
#' # --> gliomaData_0.4.tar.gz#' require(gliomaData)
#' library(gliomaData)
#' data(ge_cgh_locIGR)
#' A <- ge_cgh_locIGR$multiblocks
#' A[[3]] <- A[[3]][, -3]
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' C <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#'
#' # rgcca algorithm using the dual formulation for X1 and X2
#' # and the dual formulation for X3
#'
#' fit.rgcca <- rgcca(A,
#'   connection = C, tau = c(1, 1, 0),
#'   ncomp = c(2, 2, 1), scheme = "factorial",
#'   verbose = TRUE
#' )
#' boot.out <- bootstrap(fit.rgcca, n_boot = 50, n_cores = 2)
#' plot(boot.out, block = 1, type = "weight", ncomp = 1, n_marks = 30)
#' plot(boot.out, block = 1, type = "weight", ncomp = 2, n_marks = 30)
#' get_bootstrap(boot.out)
#'
#' # stability analysis prior bootstrap for sgcca
#' }
#' @export
#' @seealso \code{\link[RGCCA]{plot.bootstrap}},
#' \code{\link[RGCCA]{print.bootstrap}}
bootstrap <- function(rgcca_res, n_boot = 100,
                      n_cores = parallel::detectCores() - 1,
                      balanced = TRUE, keep_all_variables = FALSE,
                      verbose = TRUE) {
  if (class(rgcca_res) == "stability") {
    message("All the parameters were imported from the fitted rgcca_stability object.")
    rgcca_res <- rgcca_res$rgcca_res
  }

  if (tolower(rgcca_res$call$method) %in% c("sgcca", "spls", "spca")) {
    if (verbose) {
      message("Only selected variables were used for bootstrapping. see rgcca_stability().")
    }

    keepVar <- lapply(
      rgcca_res$a,
      function(x) unique(which(x != 0, arr.ind = TRUE)[, 1])
    )

    newBlock <- mapply(function(x, y) x[, y, drop = FALSE],
      rgcca_res$call$raw, keepVar,
      SIMPLIFY = FALSE
    )

    rgcca_res <- rgcca(newBlock,
      connection = rgcca_res$call$connection,
      superblock = rgcca_res$call$superblock,
      ncomp = rgcca_res$call$ncomp,
      bias = rgcca_res$call$bias,
      tau = 1,
      scale = rgcca_res$call$scale,
      verbose = FALSE,
      scale_block = rgcca_res$call$scale_block
    )
  }

  check_integer("n_boot", n_boot)
  check_integer("n_cores", n_cores, min = 0)

  if (n_cores == 0) n_cores <- 1

  boot_sampling <- generate_resampling(
    rgcca_res = rgcca_res,
    n_boot = n_boot,
    balanced = balanced,
    keep_all_variables = keep_all_variables,
    verbose = verbose
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
  list_res_W <- list_res_L <- list()
  for (i in 1:ndefl_max) {
    list_res_W[[i]] <- list_res_L[[i]] <- list()
    for (block in names(rgcca_res$call$blocks))
    {
      list_res_L[[i]][[block]] <- list_res_W[[i]][[block]] <-
        matrix(NA, NCOL(rgcca_res$call$blocks[[block]]), n_boot)
      rownames(list_res_L[[i]][[block]]) <-
        rownames(list_res_W[[i]][[block]]) <-
        colnames(rgcca_res$call$blocks[[block]])
    }
  }

  blocks <- NULL

  if (!verbose) {
    pbapply::pboptions(type = "none")
  } else {
    pbapply::pboptions(type = "timer")
  }

  if (Sys.info()["sysname"] == "Windows") {
    if (n_cores > 1) {
      assign("rgcca_res", rgcca_res, envir = .GlobalEnv)
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, "rgcca_res")
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
      rm("rgcca_res", envir = .GlobalEnv)
    } else {
      W <- pbapply::pblapply(
        boot_sampling$full_idx,
        function(b) {
          bootstrap_k(
            rgcca_res = rgcca_res,
            inds = b
          )
        }
      )
    }
  } else {
    W <- pbapply::pblapply(boot_sampling$full_idx,
      function(b) {
        bootstrap_k(
          rgcca_res = rgcca_res,
          inds = b
        )
      },
      cl = n_cores
    )
  }

  L <- lapply(W, `[[`, 2)
  W <- lapply(W, `[[`, 1)

  for (k in seq(n_boot)) {
    for (i in 1:ndefl_max) {
      for (j in 1:length(rgcca_res$call$blocks)) {
        block <- names(rgcca_res$call$blocks)[j]
        if (!is.null(names(W[[k]]))) {
          if (i <= NCOL(W[[k]][[block]])) {
            list_res_W[[i]][[block]][, k] <- W[[k]][[block]][, i]
            list_res_L[[i]][[block]][, k] <- L[[k]][[block]][, i]
          }
        } else {
          list_res_W[[i]][[block]][, k] <-
            rep(NA, length(list_res_W[[i]][[block]][, k]))
          list_res_L[[i]][[block]][, k] <-
            rep(NA, length(list_res_L[[i]][[block]][, k]))
          if (is.character(W[[k]]) && ((j == 1) && (i == 1))) {
            warning(paste0(
              "This bootstrap sample was removed due to zero variance variable(s): ",
              paste(W[[k]], collapse = " - ")
            ))
          }
        }
      }
    }
  }

  return(structure(list(
    bootstrap = list(W = list_res_W, L = list_res_L),
    rgcca = rgcca_res
  ),
  class = "bootstrap"
  ))
}
