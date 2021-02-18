#' Bootstrap confidence intervals and p-values
#'
#' Boostrap confidence intervals and p-values for evaluating the significancy/ 
#' stability of the block-weight vectors produce by S/RGCCA.
#' @inheritParams plot2D
#' @param n_boot Number of bootstrap iterations. Default is 100.
#' @param n_cores Number of cores for parallelization 
#' @return A list containing two objects: 'bootstrap' and 'rgcca'. 
#' 'bootstrap' is a list containing for each block, a matrix 
#' with the variables of the block in row and the block weight vector
#' calculated accross bootstrap sample in column. 'rgcca' is the original 
#' fitted rgcca object. (see  \code{\link[RGCCA]{rgcca}})
#' @examples
#' # Bootstrap confidence intervals and p-values for RGCCA 
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], 
#'               industry = Russett[, 4:5], 
#'               politic = Russett[, 6:11])
#'               
#' fit.rgcca = rgcca(blocks, ncomp= 1)
#' boot.out = bootstrap(fit.rgcca, 
#'                      n_boot = 20, n_cores = 2)
#'  
#' plot(boot.out, block = 3, comp = 1)
#' 
#' print(boot.out)
#' get_bootstrap(boot.out)
#' 
#' fit.rgcca = rgcca(blocks, type = "mcoa")
#' boot.out = bootstrap(fit.rgcca, 
#'                      n_boot = 50, n_cores = 2)
#'                      
#' plot(boot.out, block = 1)
#' 
#' \dontrun{
#' # Stability of the selected variables for SGCCA
#' # Not run: 
#' # Download the dataset's package at http://biodev.cea.fr/sgcca/.
#' # --> gliomaData_0.4.tar.gz#' require(gliomaData)
#' library(gliomaData)
#' data(ge_cgh_locIGR)
#' A <- ge_cgh_locIGR$multiblocks
#' A[[3]] = A[[3]][, -3]
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' 
#' # rgcca algorithm using the dual formulation for X1 and X2 
#' # and the dual formulation for X3
#'  
#' fit.rgcca = rgcca(A, connection = C, tau = c(1, 1, 0), 
#'                   ncomp = c(2, 2, 1), scheme = "factorial", 
#'                   verbose = TRUE)
#' boot.out = bootstrap(fit.rgcca, n_boot = 50, n_cores = 2)
#' plot(boot.out, block = 1, ncomp = 1, n_marks = 30)
#' plot(boot.out, block = 1, ncomp = 2, n_marks = 30)
#' get_bootstrap(boot.out)
#'                       
#' # sgcca algorithm
#' result.sgcca = sgcca(A, C, sparsity = c(.071,.2, 1), ncomp = 1,
#'                      scheme = "centroid", verbose = TRUE)
#'                      
#' boot.out = bootstrap(fit.sgcca, n_boot = 50, n_cores = 2)
#' }
#' @export
#' @seealso \code{\link[RGCCA]{plot.bootstrap}}, 
#' \code{\link[RGCCA]{print.bootstrap}} 
bootstrap <- function(rgcca_res, 
                      n_boot = 100, 
                      n_cores = parallel::detectCores() - 1){
    
    ndefl_max = max(rgcca_res$call$ncomp)
    list_res = list()
    for(i in 1:ndefl_max){ 
        list_res[[i]] = list()
        for(block in names(rgcca_res$blocks))
        {
            list_res[[i]][[block]] = 
              matrix(NA, dim(rgcca_res$blocks[[block]])[2], n_boot) 
            rownames( list_res[[i]][[block]]) = 
              colnames(rgcca_res$blocks[[block]])
        }
    }

    stopifnot(is(rgcca_res, "rgcca"))
    check_integer("n_boot", n_boot)
    check_integer("n_cores", n_cores, min = 0)

    if (n_cores == 0) n_cores <- 1

        blocks <- NULL

    # varlist <- c(ls(getNamespace("RGCCA")))
    # W <- RGCCA:::parallelize(
    #     varlist,
    #     seq(n_boot), 
    #     function(x) resBoot = RGCCA:::bootstrap_k(rgcca_res, type = "weight"), 
    #     n_cores = n_cores,
    #     envir = environment(),
    #     applyFunc = "parLapply",
    #     parallelization = parallelization
    #     )
    
    if(n_cores>1){
        assign("rgcca_res", rgcca_res, envir = .GlobalEnv)
        cl = parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, "rgcca_res")
        W = pbapply::pblapply(seq(n_boot), 
                function(b) bootstrap_k(rgcca_res, "weight"),
                cl = cl)
        parallel::stopCluster(cl)
        rm("rgcca_res", envir = .GlobalEnv)
    }
    else
        W = pbapply::pblapply(seq(n_boot), 
                              function(b) bootstrap_k(rgcca_res, "weight"))

    for(k in seq(n_boot)){
     for(i in 1:ndefl_max){
       for(j in 1:length(rgcca_res$blocks)){
         block = names(rgcca_res$blocks)[j]
         if(!is.null(names(W[[k]]))){
           if(i <= NCOL(W[[k]][[block]])){
             list_res[[i]][[block]][, k] = W[[k]][[block]][, i]
           }
         }
         else{
           list_res[[i]][[block]][, k] = 
               rep(NA, length(list_res[[i]][[block]][, k]))
         }
       }
     }
    }

    return(structure(list(bootstrap = list_res, 
                          rgcca = rgcca_res), 
                          class = "bootstrap"))
    

}
