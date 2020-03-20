#' Set a list of ocket for parralel package in Windows
#' @param f : a function to parralelize
#' @param nperm : a vector object for a lapply type function
#' @param varlist : character vector of names of objects to export  
#' @param envir : environment                                                  
#' @param varlist : character vector of names of objects to export                                                  
#' @param applyFunc: function to be applied
parallelize <- function(
varlist = c(),
nperm,
f,
n_cores = NULL,
envir = environment(),
applyFunc = "parSapply") {

    load_libraries("parallel")
    if (!("parallel" %in% installed.packages()[, "Package"]))
    stop("'parallel' package required and not available.")

    if(is.null(n_cores))
    n_cores <- parallel::detectCores() - 1

    if (Sys.info()["sysname"] == "Windows") {

        cl <- parallel::makeCluster(n_cores)

        parallel::clusterExport(
        cl,
        varlist,
        envir = envir
        )

        parallel::clusterEvalQ(cl, library(RGCCA))
        parallel::clusterEvalQ(cl, library(parallel))

        res <- tryCatch({
            get(applyFunc)(
            cl,
            nperm,
            f)
        }, error = function(err) stop(err$message),
        finally = {
            parallel::stopCluster(cl)
            cl <- c()
        })


    }else{

        res <- parallel::mclapply(
            nperm,
            f,
            mc.cores = n_cores)

        if (applyFunc == "parSapply")
           res <- simplify2array(res)
    }

    return(res)
}
