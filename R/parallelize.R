# Set a list of sockets for parrallel package
# @param f a function to parallelize
# @param nperm a vector object for a lapply type function
# @param varlist character vector of names of objects to export
# @param envir environment
# @param applyFunc function to be applied
# @inheritParams bootstrap
# @importFrom parallel stopCluster
# @importFrom parallel clusterExport
# @importFrom parallel clusterEvalQ
# @importFrom parallel makeCluster
# @importFrom parallel detectCores
# @importFrom parallel parLapply
# @importFrom parallel parSapply
# @importFrom parallel mclapply
#' @importFrom utils installed.packages

parallelize <- function(varlist = c(),
                        nperm,
                        f,
                        n_cores = NULL,
                        envir = environment(),
                        applyFunc = "parSapply",
                        parallelization = NULL) {
  if (is.null(parallelization)) {
    if (Sys.info()["sysname"] == "Windows") {
      parallelization <- FALSE
    } else {
      #    if( Sys.info()["sysname"] == "Windows")
      #    {
      #        message("Windows can be slow for starting parallelization.
      #                 Using parallelization=FALSE can conduct to faster
      #                 results for fast computations")
      #    }

      parallelization <- TRUE
    }
  }

  if (parallelization == TRUE && Sys.info()["sysname"] != "Windows") {
    load_libraries("parallel")
    if (!("parallel" %in% installed.packages()[, "Package"])) {
      stop_rgcca("'parallel' package required and not available.")
    }

    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }

    #  if (Sys.info()["sysname"] == "Windows") {

    #
    # cl <- parallel::makeCluster(n_cores)
    #
    # parallel::clusterExport(
    #     cl,
    #     varlist,
    #     envir = envir
    # )
    #
    #
    # parallel::clusterEvalQ(cl, library(RGCCA))
    #
    # # library(parallel)
    # parallel::clusterEvalQ(cl, library(parallel))
    #
    #      res <- tryCatch({
    #         get(applyFunc)(
    #             cl,
    #             nperm,
    #             f)
    #     }, error = function(err) stop_rgcca(err$message),
    #     finally = {
    #         parallel::stopCluster(cl)
    #         cl <- c()
    #     })
    #
    #
    # }else{
    res <- parallel::mclapply(
      nperm,
      f,
      mc.cores = n_cores
    )

    if (applyFunc == "parSapply") {
      res <- simplify2array(res)
    }

    #  }
  }
  if (Sys.info()["sysname"] == "Windows" || parallelization == FALSE) {
    res <- lapply(
      nperm,
      f
    )
    if (applyFunc == "parSapply") {
      res <- simplify2array(res)
    }
  }
  if (parallelization == "for") {
    res <- NULL
    for (i in 1:length(nperm))
    {
      res[[i]] <- f(nperm[i])
    }
    if (applyFunc == "parSapply") {
      res <- simplify2array(res)
    }
  }

  return(res)
}
