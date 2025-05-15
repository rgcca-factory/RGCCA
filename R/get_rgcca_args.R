#' Get the rgcca arguments from a fitted object and default arguments.
#' Modify arguments based on the provided configuration if needed.
#' @noRd
get_rgcca_args <- function(object, default_args = list()) {
  if (any(class(object) %in% c("rgcca", "rgcca_permutation", "rgcca_cv"))) {
    opt <- object$opt
    rgcca_args <- object$call

    if (any(class(object) %in% c("rgcca_permutation", "rgcca_cv"))) {
      if (object$par_type == "tau") {
        rgcca_args$tau <- object$best_params
      }
      if (object$par_type == "ncomp") {
        rgcca_args$ncomp <- object$best_params
      }
      if (object$par_type == "sparsity") {
        rgcca_args$sparsity <- object$best_params
      }
    }
  } else {
    rgcca_args <- list(
      tau = default_args$tau,
      tol = default_args$tol,
      init = tolower(default_args$init),
      bias = default_args$bias,
      quiet = default_args$quiet,
      scale = default_args$scale,
      ncomp = default_args$ncomp,
      blocks = default_args$blocks,
      scheme = default_args$scheme,
      method = tolower(default_args$method),
      verbose = default_args$verbose,
      sparsity = default_args$sparsity,
      response = default_args$response,
      NA_method = tolower(default_args$NA_method),
      comp_orth = default_args$comp_orth,
      n_iter_max = default_args$n_iter_max,
      connection = default_args$connection,
      superblock = default_args$superblock,
      scale_block = default_args$scale_block,
      lambda = default_args$lambda,
      graph_laplacians = default_args$graph_laplacians
    )

    rgcca_args$init <- check_char(rgcca_args$init, "init", c("svd", "random"))
    rgcca_args$NA_method <- check_char(
      rgcca_args$NA_method, "NA_method", c("na.ignore", "na.omit")
    )
    if (!is.logical(rgcca_args$scale_block)) {
      rgcca_args$scale_block <- tolower(rgcca_args$scale_block)
      rgcca_args$scale_block <- check_char(
        rgcca_args$scale_block, "scale_block", c("inertia", "lambda1")
      )
    }

    rgcca_args$blocks <- check_blocks(
      rgcca_args$blocks, add_NAlines = TRUE,
      quiet = rgcca_args$quiet, response = rgcca_args$response
    )

    check_integer("tol", rgcca_args$tol, float = TRUE, min = 0)
    check_integer("n_iter_max", rgcca_args$n_iter_max, min = 1)
    for (i in c(
      "superblock", "verbose", "scale", "bias", "quiet", "comp_orth"
    )) {
      check_boolean(i, rgcca_args[[i]])
    }
    rgcca_args$graph_laplacians <- check_laplacians(rgcca_args$graph_laplacians,
      rgcca_args$blocks)

    rgcca_args$tau <- elongate_arg(rgcca_args$tau, rgcca_args$blocks)
    rgcca_args$ncomp <- elongate_arg(rgcca_args$ncomp, rgcca_args$blocks)
    rgcca_args$sparsity <- elongate_arg(rgcca_args$sparsity, rgcca_args$blocks)
    # a lambda that is diff. from 0 does not matter if lap. is NULL
    rgcca_args$lambda <- elongate_arg(rgcca_args$lambda, rgcca_args$blocks)

    ### Get last parameters based on the method
    tmp <- select_analysis(rgcca_args, rgcca_args$blocks)
    opt <- tmp$opt
    rgcca_args <- tmp$rgcca_args

    # Change penalty to 0 if there is a univariate disjunctive block response
    opt$disjunction <- !is.null(rgcca_args$response) &&
      is.character(rgcca_args$blocks[[rgcca_args$response]])

    if (opt$disjunction) {
      if (is.matrix(rgcca_args[[opt$param]])) {
        rgcca_args[[opt$param]][, rgcca_args$response] <- 0
      } else {
        rgcca_args[[opt$param]][rgcca_args$response] <- 0
      }
    }
  }
  rgcca_args$quiet <- TRUE
  rgcca_args$verbose <- FALSE
  return(list(opt = opt, rgcca_args = rgcca_args))
}
