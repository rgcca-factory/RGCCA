#' Bootstrap confidence intervals and p-values
#'
#' Boostrap confidence intervals and p-values for evaluating the significancy/
#' stability of the block-weight vectors produce by S/RGCCA.
#' @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})
#' @param n_boot Number of bootstrap samples. Default: 100.
#' @param n_cores Number of cores for parallelization.
#' @return A list containing two objects: 'bootstrap' and 'rgcca'.
#' 'bootstrap' is a list containing for each block, a matrix
#' with the variables of the block in row and the block weight vector
#' calculated accross bootstrap sample in column. 'rgcca' is the fitted rgcca
#' object obtained from the original data. (see  \code{\link[RGCCA]{rgcca}})
#' @examples
#' # Bootstrap confidence intervals and p-values for RGCCA
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)],
#'               industry = Russett[, 4:5],
#'               politic = Russett[, 6:11])
#'
#' fit.rgcca = rgcca(blocks, ncomp= c(2, 1, 2))
#' boot.out = bootstrap(fit.rgcca, n_boot = 20, n_cores = 2)
#'
#' plot(boot.out, type = "weight", block = 3, comp = 1)
#'
#' print(boot.out, comp =2)
#'  get_bootstrap(boot.out, block = 1, comp = 1)
#'
#' fit.rgcca = rgcca(blocks, method = "mcoa")
#' boot.out = bootstrap(fit.rgcca, n_boot = 50, n_cores = 2)
#'
#' plot(boot.out, type = "weight", block = 4)
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
                      n_cores = parallel::detectCores() - 1){

    if(class(rgcca_res)=="stability")
    {
        message("All the parameters were imported from the fitted rgcca_stability object.")
        rgcca_res = rgcca_res$rgcca_res
    }
    
    stopifnot(is(rgcca_res, "rgcca"))
    if(tolower(rgcca_res$call$method)%in%c("sgcca", "spls", "spca"))
        stop_rgcca("for sparse models, rgcca_stability() applies before bootstrapping")
    check_integer("n_boot", n_boot)
    check_integer("n_cores", n_cores, min = 0)
    
    if (n_cores == 0) n_cores <- 1

    ndefl_max = max(rgcca_res$call$ncomp)
    list_res_W = list_res_L = list()
    for(i in 1:ndefl_max){
        list_res_W[[i]] = list_res_L[[i]] = list()
        for(block in names(rgcca_res$call$blocks))
        {
            list_res_L[[i]][[block]] = list_res_W[[i]][[block]] =
                matrix(NA, NCOL(rgcca_res$call$blocks[[block]]), n_boot)
            rownames(list_res_L[[i]][[block]]) =
                rownames(list_res_W[[i]][[block]]) =
                colnames(rgcca_res$call$blocks[[block]])
        }
    }

    blocks <- NULL

    if( Sys.info()["sysname"] == "Windows"){
    if(n_cores>1){
        assign("rgcca_res", rgcca_res, envir = .GlobalEnv)
        cl = parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, "rgcca_res")
        W = pbapply::pblapply(seq(n_boot),
                function(b) bootstrap_k(rgcca_res),
                cl = cl)
        parallel::stopCluster(cl)
        rm("rgcca_res", envir = .GlobalEnv)
    }
    else
        W = pbapply::pblapply(seq(n_boot),
                              function(b) bootstrap_k(rgcca_res))
    }else{
        W = pbapply::pblapply(seq(n_boot),
                              function(b) bootstrap_k(rgcca_res),
                              cl = n_cores)
    }
    
    if (sum(unlist(lapply(W, is.na))) != 0){
        boot_to_rm    = which(!unlist(lapply(W, is.list)))
        W[boot_to_rm] = NULL
        list_res      = lapply(list_res, function(x) lapply(x, function(y) y[, -boot_to_rm]))
        warning("Even after multiple resampling, some bootstrap samples could
                 not overcome the projection issue mentioned in the previous 
                 warnings. Hence, the number of bootstrap samples is decreased 
                 from ", n_boot, " to ", length(W), ".")
        n_boot        = length(W)
    }

    L = lapply(W, `[[`, 2)
    W = lapply(W, `[[`, 1)

    for(k in seq(n_boot)){
     for(i in 1:ndefl_max){
       for(j in 1:length(rgcca_res$call$blocks)){
         block = names(rgcca_res$call$blocks)[j]
         if(!is.null(names(W[[k]]))){
           if(i <= NCOL(W[[k]][[block]])){
             list_res_W[[i]][[block]][, k] = W[[k]][[block]][, i]
             list_res_L[[i]][[block]][, k] = L[[k]][[block]][, i]
           }
         }
         else{
           list_res_W[[i]][[block]][, k] =
               rep(NA, length(list_res_W[[i]][[block]][, k]))
           list_res_L[[i]][[block]][, k] =
               rep(NA, length(list_res_L[[i]][[block]][, k]))
         }
       }
     }
    }

    return(structure(list(bootstrap = list(W = list_res_W, L = list_res_L),
                          rgcca = rgcca_res),
                          class = "bootstrap"))

}
