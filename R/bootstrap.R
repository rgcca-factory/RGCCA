#' Bootstrap confidence intervals and p-values
#'
#' Boostrap confidence intervals and p-values for evaluating the significancy/ 
#' stability of the block-weight vectors produce by S/RGCCA.
#' @inheritParams plot2D
#' @param n_boot Number of bootstrap iterations. Default is 100.
#' @param n_cores Number of cores for parallelization 
#' @param parallelization if TRUE, the bootstrap is processed in parallel. If 
#' parallelization = NULL (default), parallelization is always performed except 
#' for Windows if length(n_boot) < 10.
#' @return A list containing two objects: 'bootstrap' and 'rgcca'. 
#' 'bootstrap' is a list containing for each block, a matrix 
#' with the variables of the block in row and the block weight vector
#' calculated accross bootstrap sample in column. 'rgcca' is the original 
#' fitted rgcca object. (see  \code{\link[RGCCA]{rgcca}})
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], 
#'               industry = Russett[, 4:5], 
#'               politic = Russett[, c(6:9, 11)]
#'               )
#'               
#' fit.rgcca = rgcca(blocks)
#' boot_rgcca = bootstrap(fit.rgcca, n_boot = 50)
#' plot(boot_rgcca, block = 3)
#' 
#' @export
#' @seealso \code{\link[RGCCA]{plot.bootstrap}} , 
#' \code{\link[RGCCA]{print.bootstrap}} 
bootstrap <- function(
    rgcca_res,
    n_boot = 100,
    n_cores = parallel::detectCores() - 1,
    parallelization = NULL) {
    ndefl_max = max(rgcca_res$call$ncomp)
    list_res = list()
    for(i in 1:ndefl_max)
    { 
        list_res[[i]] = list()
        for(block in names(rgcca_res$call$blocks))
        {
            list_res[[i]][[block]] = 
              matrix(NA, dim(rgcca_res$call$blocks[[block]])[2], n_boot) 
            rownames( list_res[[i]][[block]]) = 
              colnames(rgcca_res$call$blocks[[block]])
        }
    }

    stopifnot(is(rgcca_res, "rgcca"))
    if (!is.null(parallelization))
        check_boolean("parallelization", parallelization)
    check_integer("n_boot", n_boot)
    check_integer("n_cores", n_cores, min = 0)

    if (n_cores == 0)
        n_cores <- 1

    message("Bootstrap in progress...", appendLF = F)

    blocks <- NULL

    varlist <- c(ls(getNamespace("RGCCA")))
    # get the parameter dot-dot-dot
    # args_values <- list(...)
    # args_names <- names(args_values)
    # n <- args_values
    # if (!is.null(n))
    #     n <- seq(length(args_values))
    # for (i in n) {
    #     if (!is.null(args_names[i])) {
    #         # dynamically asssign these values
    #         assign(args_names[i], args_values[[i]])
    #         # send them to the clusters to parallelize
    #         varlist <- c(varlist, args_names[i])
    #         # without this procedure rgcca_cv_k(rgcca_res, blocks = blocks2)
    #         # or rgcca_cv_k(rgcca_res, blocks = lapply(blocks, scale)
    #         # does not work.
    #     }
    # }
 
    W <- parallelize(
        varlist,
        seq(n_boot), 
        function(x) {
            resBoot = bootstrap_k(rgcca_res, type = "weight")
        }
        , 
        n_cores = n_cores,
        envir = environment(),
        applyFunc = "parLapply",
        parallelization = parallelization
        )

       for(k in seq(n_boot))
       {
               for(i in 1:ndefl_max)
              {
                   for(j in 1:length(rgcca_res$call$blocks))
                   {
                       block=names(rgcca_res$call$blocks)[j]
                   
                       if(!is.null(names(W[[k]])))
                       {
                           if(i<=dim(W[[k]][[block]])[2])
                          {
                             list_res[[i]][[block]][, k] = W[[k]][[block]][, i]
                
                          } 
                       }
                       else
                       {
                           list_res[[i]][[block]][, k] = 
                             rep(NA, length(list_res[[i]][[block]][, k]))
                       }
                   }
                }
           }
      
       
      message("OK.")
      
      return(structure(list(bootstrap = list_res, 
                            rgcca = rgcca_res), 
                       class = "bootstrap")
             )
}
