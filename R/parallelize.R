#' Set a list of sockets for parralel package
#' @param f a function to parralelize
#' @param nperm a vector object for a lapply type function
#' @param varlist character vector of names of objects to export  
#' @param envir environment                                               
#' @param applyFunc function to be applied
#' @inheritParams bootstrap
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel parLapply
#' @importFrom parallel parSapply
parallelize <- function(
    varlist = c(),
    nperm,
    f,
    n_cores = NULL,
    envir = environment(),
    applyFunc = "parSapply") {

#    load_libraries("parallel")
#   if (!("parallel" %in% installed.packages()[, "Package"]))
#    stop_rgcca("'parallel' package required and not available.")

    if(is.null(n_cores))
    n_cores <- detectCores() - 1

    if (Sys.info()["sysname"] == "Windows") {

      
            cl <- makeCluster(n_cores)
            
            clusterExport(
                cl,
                varlist,
                envir = envir
            )
            
            
            clusterEvalQ(cl, library(RGCCA))
      
            # library(parallel)
            clusterEvalQ(cl, library(parallel))
            # print("all okay")
            res <- tryCatch({
                get(applyFunc)(
                    cl,
                    nperm,
                    f)
            }, error = function(err) stop_rgcca(err$message),
            finally = {
                stopCluster(cl)
                cl <- c()
            })
            
      
        

    }else{

        res <- mclapply(
            nperm,
            f,
            mc.cores = n_cores)

        if (applyFunc == "parSapply")
           res <- simplify2array(res)
    }

    return(res)
}
