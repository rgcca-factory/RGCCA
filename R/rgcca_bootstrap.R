#' Bootstrap confidence intervals and p-values
#'
#' Bootstrap confidence intervals and p-values for evaluating the significance/
#' stability of the block-weight vectors produce by S/RGCCA.
#' @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})
#' @param n_boot Number of bootstrap samples. Default: 100.
#' @param n_cores Number of cores for parallelization.
#' @param balanced A boolean indicating if a balanced bootstrap procedure is
#' performed or not (default is TRUE).
#' @param keep_all_variables A boolean indicating if all variables have to be
#' kept even when some of them have null variance for at least one bootstrap
#' sample (default is FALSE).
#' @param verbose Logical value indicating if the progress of the
#' bootstrap procedure is reported.
#' @return A list containing two objects: 'bootstrap' and 'rgcca'.
#' 'bootstrap' is a list containing for each block, a matrix
#' with the variables of the block in row and the block weight vector
#' calculated across bootstrap sample in column. 'rgcca' is the fitted rgcca
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
#' boot.out <- rgcca_bootstrap(fit.rgcca, n_boot = 20, n_cores = 2)
#'
#' plot(boot.out, type = "weight", block = 3, comp = 1)
#'
#' print(boot.out, comp = 2, block = 1)
#'
#' fit.rgcca <- rgcca(blocks, method = "mcoa")
#' boot.out <- rgcca_bootstrap(fit.rgcca, n_boot = 50, n_cores = 2)
#'
#' plot(boot.out, type = "weight", block = 1)
#' \dontrun{
#' # Stability of the selected variables for SGCCA
#' # Download the dataset's package at http://biodev.cea.fr/sgcca/ and install
#' # it from the package archive file.
#' # You can do it with the following R commands:
#' if (!("gliomaData" %in% rownames(installed.packages()))) {
#'   destfile <- tempfile()
#'   download.file("http://biodev.cea.fr/sgcca/gliomaData_0.4.tar.gz", destfile)
#'   install.packages(destfile, repos = NULL, type = "source")
#' }
#'
#' data("ge_cgh_locIGR", package = "gliomaData")
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
#' boot.out <- rgcca_bootstrap(fit.rgcca, n_boot = 50, n_cores = 2)
#' plot(boot.out, block = 1, type = "weight", ncomp = 1, n_marks = 30)
#' plot(boot.out, block = 1, type = "weight", ncomp = 2, n_marks = 30)
#'
#' # stability analysis prior bootstrap for sgcca
#' }
#' @export
#' @seealso \code{\link[RGCCA]{plot.bootstrap}},
#' \code{\link[RGCCA]{print.bootstrap}}
rgcca_bootstrap <- function(rgcca_res, n_boot = 100,
                            n_cores = 1,
                            balanced = TRUE, keep_all_variables = FALSE,
                            verbose = TRUE) {
  stability <- is(rgcca_res, "stability")
  if (stability) {
    message(
      "All the parameters were imported from the fitted rgcca_stability",
      " object."
    )
    rgcca_res <- rgcca_res$rgcca_res
  }

  # If sparse model, we perform bootstrap only on the selected variables
  if (!stability && tolower(rgcca_res$call$method) %in% sparse_methods()) {
    if (verbose) {
      message(
        "Only selected variables were used for bootstrapping. see ",
        "rgcca_stability()."
      )
    }

    # Remove superblock variables from keep_var as the superblock is generated
    # from the kept variables
    J <- length(rgcca_res$call$blocks)
    keep_var <- lapply(
      rgcca_res$a[-(J + 1)],
      function(x) unique(which(x != 0, arr.ind = TRUE)[, 1])
    )
    if (rgcca_res$opt$disjunction) {
      keep_var[[rgcca_res$call$response]] <- 1
    }

    rgcca_res$call$blocks <- Map(
      function(x, y) x[, y, drop = FALSE], rgcca_res$call$blocks, keep_var
    )
    rgcca_res$call$tau <-
      rgcca_res$call$sparsity <- rep(1, length(rgcca_res$blocks))

    rgcca_res <- rgcca(rgcca_res)
  }

  check_integer("n_boot", n_boot)

  boot_sampling <- generate_resampling(
    rgcca_res = rgcca_res,
    n_boot = n_boot,
    balanced = balanced,
    keep_all_variables = keep_all_variables,
    verbose = verbose
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
        inds = b
      )
    },
    n_cores = n_cores, verbose = verbose
  )

  res <- format_bootstrap_list(W, rgcca_res)

  if (verbose) {
    message("Computation of bootstrap statistics...")
  }
  stats <- rgcca_bootstrap_stats(res, rgcca_res, verbose = verbose)

  return(structure(list(
    n_boot = n_boot, rgcca = rgcca_res,
    bootstrap = res, stats = data.frame(stats)
  ),
  class = "bootstrap"
  ))
}
