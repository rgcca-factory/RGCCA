#' Set a list of ocket for parralel package in Windows
#' f : a function to parralelize
#' nperm : a vector object for a lapply type function
#' varlist : character vector of names of objects to export                                                  
parallelize <- function(
    varlist,
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
        }, error = function(e) print(e), 
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
