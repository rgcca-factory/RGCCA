#' Set a list of sockets for parralel package
#' @param f a function to parralelize
#' @param nperm a vector object for a lapply type function
#' @param varlist character vector of names of objects to export  
#' @param envir environment                                               
#' @param applyFunc function to be applied
#' @param para if TRUE parallelization is run, if FALSE, no parallelisation is run. If NULL (default) parallelization is always used except for Windows in case of length(nperm)<10
#' @inheritParams bootstrap
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel parLapply
#' @importFrom parallel parSapply
#' @importFrom parallel mclapply
parallelize <- function(
    varlist = c(),
    nperm,
    f,
    n_cores = NULL,
    envir = environment(),
    applyFunc = "parSapply", 
    para=NULL) {
    if(is.null(para))
    {
        if( Sys.info()["sysname"] == "Windows"& length(nperm)<10)
        {
            para=FALSE
        }
        else
        {
            message("Windows can be slow for starting parallelization. Using para=FALSE can conduct to faster results for light calculations")
            para=TRUE
        }
       
    }
    
    if(para)
    {
        load_libraries("parallel")
        if (!("parallel" %in% installed.packages()[, "Package"]))
            stop_rgcca("'parallel' package required and not available.")
        
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
            
            # library(parallel)
            parallel::clusterEvalQ(cl, library(parallel))
            # print("all okay")
            res <- tryCatch({
                get(applyFunc)(
                    cl,
                    nperm,
                    f)
            }, error = function(err) stop_rgcca(err$message),
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
        
    }
    else
    {
        res=lapply(
            nperm,
            f)
        if (applyFunc == "parSapply")
            res <- simplify2array(res)
    }
   
    return(res)
}
