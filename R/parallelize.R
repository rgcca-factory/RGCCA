#' Set a list of ocket for parralel package in Windows
#' f : a function to parralelize
#' nperm : a vector object for a lapply type function
#' varlist : character vector of names of objects to export                                                                                                                                                                                                                                                                                                                                                                          
parallelize <- function(varlist, nperm, f, n_cores = parallel::detectCores() - 1, envir = environment()){

    if (Sys.info()["sysname"] == "Windows") {

        cl <- parallel::makeCluster(n_cores)
        
        parallel::clusterExport(
            cl,
            varlist,
            envir = envir
        )
 
        parallel::clusterEvalQ(cl, library(RGCCA))

        res <- tryCatch({
            parallel::parSapply(
                cl,
                nperm,
                f)
        }, error = function(e) print(e), 
        finally = {
            parallel::stopCluster(cl)
            cl <- c()
        })


    }else{

        res <- simplify2array(
            parallel::mclapply(
                nperm,
                f,
                mc.cores = n_cores))
    }

    return(res)
}
