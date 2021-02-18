#' Tune the S/RGCCA hyper-parameters by permutation.
#' 
#' This function can be used to automatically select the hyper-parameters 
#' (amount of sparsity for sgcca or shrinkage parameters for RGCCA)
#' A permutation based strategy very similar to the one proposed in 
#' (Witten et al, 2009) is proposed.
#' 
#' @details 
#' The tuning parameters are selected using the permutation scheme proposed in 
#' (Witten et al, 2009). For each candidate tuning parameter value, the 
#' following is performed: 
#' 
#' (1) Repeat the following n_perms times (for n_perms large): \cr
#'    \verb{    }(a) The samples in \eqn{X_1},..., \eqn{X_J} are randomly 
#'    permuted blocks: \eqn{X_1^*},..., \eqn{X_J^*}. \cr
#'    \verb{    }(b) S/RGCCA is run on the permuted data sets \eqn{X_1^*},..., 
#'       \eqn{X_J^*} to get canonical variates \eqn{a_1^*},..., \eqn{a_J^*}.\cr
#'    \verb{    }(c) Record t* = sum_(j,k) c_jk g(Cov(X_j^*a_j^*, X_k^*a_k^*). 
#'       
#' (2) Sparse CCA is run on the blocks \eqn{X_1},..., \eqn{X_J} to obtain 
#'     canonical variates \eqn{a_1},..., \eqn{a_J}. 
#' 
#' (3) Record t = sum_(j,k) c_jk g(Cov(X_ja_j, X_ka_k). 
#' 
#' (4) The resulting p-value is given by $mean(t* > t)$; that is, the fraction 
#' of t* that exceed the value of t obtained from the real data. 
#' 
#' Then, choose the tuning parameter value that gives the smallest value in 
#' Step 4.
#' 
#' This function only selects tuning parameters for the first deflation stage 
#' of S/RGCCA. By default, this function performs a one-dimensional 
#' search in tuning parameter space.
#' 
#' @inheritParams set_connection
#' @inheritParams bootstrap
#' @inheritParams rgcca
#' @inheritParams plot2D
#' @param par_type A character string indicating the parameter to tune between 
#' "sparsity" and "tau"
#' @param par_length An integer indicating the number of sets of parameters to 
#' be tested (if par_value = NULL). The parameters are uniformly distributed.
#' @param par_value It could be either (i) A matrix of dimension IxJ (where 
#' J the number of blocks and I the number of combinations to be tested), or 
#' (ii) a vector of length J length ) or 
#' a numeric value giving sets of penalties (tau for RGCCA, sparsity for SGCCA) 
#' to be tested, one row by combination. By default, it takes 10 sets between 
#' min values (0 for RGCCA and 1/sqrt(ncol(Xj)) for SGCCA) and 1.
#' @param n_perms Number of permutations for each set of constraints (default 
#' is 20).
#' @param parallelization logical value. If TRUE (default value), the 
#' permutation procedure is parallelized
#' @return \item{zstat}{The vector of Z-statistics, one per set of tuning 
#' parameters}
#' @return \item{bestpenalties}{The set of tuning parameters that gives the 
#' highest Z-statistics}
#' @return \item{permcrit}{Matrix of permuted S/RGCCA criteria. The ith row of 
#' permcrit contains the n_perms values of the permuted S/RGCCA criteria 
#' obtained for each set of tuning parameters.}
#' @return \item{means}{A vector that contains, for each set of tuning 
#' parameters, the mean of the permuted R/SGCCA criteria}
#' @return \item{sds}{A vector that contains, for each set of tuning 
#' parameters, the standard deviation of the permuted R/SGCCA criteria}
#' @return \item{crit}{A vector that contains, for each set of tuning 
#' parameters, the value of the R/SGCCA criteria obtained from the original 
#' blocks}
#' @return \item{pval}{The vector of p-values, one per set of tuning parameters}
#' @return \item{penalties}{A matrix giving, the set of tuning paramaters 
#' considered during the permutation procedure (tau or sparsity).}
#' @references Witten, D. M., Tibshirani, R., & Hastie, T. (2009). A penalized 
#' matrix decomposition, with applications to sparse principal components and 
#' canonical correlation analysis. Biostatistics, 10(3), 515-534.
#' @examples
#' ####################################
#' # Permutation based strategy for   #
#' # determining the best shrinkage   # 
#' # parameters (par_type = "tau")    #
#' ####################################
#' 
#' data(Russett)
#' blocks = list(agriculture = Russett[, seq(3)], 
#'               industry = Russett[, 4:5], 
#'               politic = Russett[, 6:11] )
#'               
#' C <-  matrix(c(0, 0, 1, 
#'                0, 0, 1, 
#'                1, 1, 0), 3, 3)
#' 
#' # default value: 10 vectors from rep(0, length(blocks))
#' # to rep(1, length(blocks)), uniformly distributed.
#' 
#' fit = rgcca_permutation(blocks, connection = C,
#'                         par_type = "tau", 
#'                         par_length = 10, n_perms = 20,
#'                         n_cores = 1)#parallel::detectCores() - 1
#'                         
#' print(fit)
#' plot(fit)
#' fit$bestpenalties
#' 
#' # It is possible to define explicitly K combinations of shrinkage
#' # parameters to be tested and in that case a matrix of dimension KxJ is 
#' # required. Each row of this matrix corresponds to one specific set of 
#' # shrinkage parameters. 
#' par_value = matrix(c(0, 0, 0, 
#'                       1, 1, 0,
#'                       0.5, 0.5, 0.5,
#'                       sapply(blocks, RGCCA:::tau.estimate), 
#'                       1, 1, 1), 5, 3, byrow = TRUE)
#' 
#'                                                                  
#' perm.out <- rgcca_permutation(blocks, connection = C,
#'                          par_type = "tau", 
#'                          par_value = par_value, 
#'                          n_perms = 5, n_cores = 1)
#' 
#' print(perm.out)
#' plot(perm.out)
#' 
#' # with superblock
#'
#' perm.out <- rgcca_permutation(blocks, par_type = "tau", 
#'                                  superblock = TRUE,
#'                                  scale = TRUE, scale_block = FALSE,
#'                                  n_perms = 5, n_cores = 1)
#' 
#' print(perm.out)
#' plot(perm.out)
#' 
#' # used a fitted rgcca_permutation object as input of the rgcca function
#' fit.rgcca <- rgcca(perm.out)
#' fit.rgcca$call$tau
#' fit.rgcca$call$scale
#' fit.rgcca$call$scale_block
#' 
#' ######################################
#' # Permutation based strategy for     #
#' # determining the best sparsity      # 
#' # parameters (par_type = "sparsity") #
#' ######################################
#' 
#' # defaut value: 10 vectors from minimum values 
#' # (1/sqrt(ncol(X1)), ..., 1/sqrt(ncol(XJ)) 
#' # to rep(1, J), uniformly distributed.
#' 
#' perm.out <- rgcca_permutation(blocks, par_type = "sparsity", 
#'                          n_perms = 50, n_cores = 1)
#'                          
#' print(perm.out)
#' plot(perm.out)
#' perm.out$bestpenalties  
#' 
#' # when par_value is a vector of length J. Each element of the vector 
#' # indicates the maximum value of sparsity to be considered for each block.
#' # par_length (default value = 10) vectors from minimum values 
#' # (1/sqrt(ncol(X1)), ..., 1/sqrt(ncol(XJ)) to maximum values, uniformly 
#' # distributed, are then considered. 
#' 
#' perm.out <- rgcca_permutation(blocks, connection = C, 
#'                          par_type = "sparsity", 
#'                          par_value = c(0.6, 0.75, 0.5), 
#'                          par_length = 7, n_perms = 20, 
#'                          n_cores = 1, tol = 1e-3)
#'                          
#' 	print(perm.out)
#' 	plot(perm.out)
#' 	perm.out$bestpenalties
#' 	
#' # when par_value is a scalar, the same maximum value is applied 
#' # for each block
#' 
#' perm.outt <- rgcca_permutation(blocks, connection = C, 
#'                          par_type = "sparsity", 
#'                          par_value = 0.8, par_length = 5, 
#'                          n_perms = 10, n_cores = 1)
#'    
#' perm.out$penalties    
#' \dontrun{
#' ######################################
#' # speed up the permutation procedure #
#' ######################################
#' 
#' # The rgcca_permutation function can be quite time-consuming. Since  
#' # approximate estimates of the block weight vectors are acceptable in this 
#' # case, it is possible to reduce the value of the tolerance (tol argument) 
#' # of the RGCCA algorithm to speed up the permutation procedure.
#' # 
#' require(gliomaData)
#' data(ge_cgh_locIGR)
#' A <- ge_cgh_locIGR$multiblocks
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' A[[3]] = A[[3]][, -3]
#' C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' 
#' # check dimensions of the blocks
#' sapply(A, dim)
#' 
#' par_value = matrix(c(seq(0.1, 1, by = 0.1),
#'                      seq(0.1, 1, by = 0.1),
#'                      rep(0, 10)), 10, 3, byrow = FALSE)
#'                        
#' fit <- rgcca_permutation(A, connection = C, 
#'                           par_type = "tau", 
#'                           par_value = par_value,
#'                           par_length = 10,
#'                           n_perms = 10, n_cores = 1, tol = 1e-2)
#'                           )
#' print(fit)
#' plot(fit)
#' }
#' 
#' @export
rgcca_permutation <- function(blocks, par_type, par_value = NULL,
                              par_length = 10, n_perms = 20, 
                              n_cores = parallel::detectCores() - 1,
                              quiet = TRUE, scale = TRUE, scale_block = TRUE, 
                              type=NULL, connection = 1 - diag(length(blocks)),
                              scheme = "factorial", 
                              ncomp = rep(1, length(blocks)), 
                              tau = rep(1, length(blocks)), 
                              sparsity = rep(1, length(blocks)), 
                              init = "svd", bias = TRUE, tol = 1e-8, 
                              response = NULL, superblock = FALSE, 
                              method = "nipals", rgcca_res = NULL, 
                              parallelization = NULL){

    if (!is.null(rgcca_res)) {
        stopifnot(is(rgcca_res, "rgcca"))
        message("All parameters were imported from a fitted rgcca object.")
        scale_block <- rgcca_res$call$scale_block
        scale <- rgcca_res$call$scale
        scheme <- rgcca_res$call$scheme
        response <- rgcca_res$call$response
        tol <- rgcca_res$call$tol
        method <- rgcca_res$call$method
        init <- rgcca_res$call$init
        bias <- rgcca_res$call$bias
        blocks <- rgcca_res$call$raw
        superblock <- rgcca_res$call$superblock
        connection <- rgcca_res$call$connection
        tau <- rgcca_res$call$tau
        ncomp <- rgcca_res$call$ncomp
        sparsity <- rgcca_res$call$sparsity
        type <- rgcca_res$call$type
        superblock <- rgcca_res$call$superblock
    }
   
    check_integer("n_perms", n_perms)
    check_integer("par_length", n_perms)
    check_integer("par_value", n_perms, min = 0)
    check_integer("n_cores", n_cores, min = 0)
    match.arg(par_type, c("tau", "sparsity"))
    min_spars <- NULL
    
    if (par_type == "sparsity") type2 <- "sgcca"
    if(par_type == "tau") type2 <- "rgcca"
    if(is.null(type)){type=type2}
    if (length(blocks) == 1) 
      stop_rgcca("Permutation required more than one block.\n")

    if(!superblock){ 
      ncols <- sapply(blocks, NCOL)
      call=list(type = type, par_type = par_type, par_value = par_value, 
                n_perms = n_perms, quiet = quiet, connection = connection, 
                method=method, tol=tol, scheme = scheme, 
                scale = scale, scale_block = scale_block, 
                superblock = superblock, blocks=blocks)
    }else{
      ncol_block = sapply(blocks, NCOL)
      ncols <- c(ncol_block, sum(ncol_block))
      names(ncols)=c(names(ncol_block), "superblock")
      J=length(ncols)
      matConnection=matrix(0,J,J);
      matConnection[1:(J-1),J]=1;matConnection[J,1:(J-1)]=1
      rownames(matConnection) = colnames(matConnection) = names(ncols)
      call=list(type = type, par_type = par_type, n_perms = n_perms, 
                quiet = quiet, connection = matConnection, 
                method=method, tol=tol, scheme = scheme, 
                scale = scale, scale_block = scale_block, 
                superblock = superblock, blocks=blocks)
    }
 
    set_spars <- function(max = 1){
      if (length(max) == 1) f <- quote(max)
      else f <- quote(max[x])
      sapply(seq(min_spars), 
        function(x) seq(eval(f), min_spars[x], len = par_length))
    }
    
    set_penalty <- function(){
      # Selecting the minimal value
      if (par_type == "sparsity"){
        if (!tolower(type) %in% c("spls", "sgcca"))
          warning("par_type = 'sparsity' but sparsity is not required... SGCCA is performed.")
          type <<- "sgcca"
          min_spars <<- sapply(ncols, function(x) 1 / sqrt(x))
      }
      else{
        if (tolower(type) %in% c("spls", "sgcca"))
          warning("par_type = 'tau' but sparsity is required... RGCCA is performed.")
          type <<- "rgcca"
          min_spars <<- sapply(ncols, function(x) 0)
      }
      
      
      if (is.null(par_value)) par_value <- set_spars()
      else if("data.frame"%in%class(par_value) || "matrix"%in%class(par_value)) 
        {
          if(par_type=="tau"){
            par_value <- t(sapply(seq(NROW(par_value)), 
                                  function(x) check_tau(par_value[x, ], 
                                                        blocks, 
                                                        type = type, 
                                                        superblock = superblock)
                                  )
                          )
          }
                     
        }
        else{#when par_value is a vector
          if (any(par_value < min_spars))
            stop_rgcca(paste0("par_value should be upper than : ", 
                                  paste0(round(min_spars, 2), collapse = ",")))
          if(par_type == "tau"){
            par_value <- check_tau(par_value, blocks, type = type, 
                                     superblock = superblock)
            par_value <- set_spars(max = par_value)
          }
          if(par_type == "sparsity"){par_value <- set_spars(max = par_value)}
        }

        if(superblock){coln = c(names(blocks), "superblock")}
        else{coln = names(blocks)}
        if(is.null(dim(par_value))){par_value = matrix(par_value, nrow = 1)}
        colnames(par_value) <- coln
        return(list(par_type, par_value))
    }

    switch(par_type, 
           "sparsity" = par <- set_penalty(), 
           "tau" = par <- set_penalty())
    
    par_value_parallel = 
      matrix(apply(par[[2]], 1, function(x) rep(x, n_perms + 1)), 
             ncol = NCOL(par[[2]]), byrow = TRUE)
    
    perm_parallel = rep(rep(c(FALSE, TRUE), c(1, n_perms)), NROW(par[[2]]))
    
    if(n_cores>1){
      call_perm = list(blocks, par[[1]], par_value_parallel, perm_parallel,
                     type, quiet, superblock, scheme,
                     tol, scale, scale_block, connection, method, 
                     response, bias, init, ncomp, tau, sparsity)
      assign("call_perm", call_perm, envir = .GlobalEnv)
      cl = parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, "call_perm")
      W = pbapply::pbsapply(seq(length(call_perm[[4]])), 
                          function(b) rgcca_permutation_k(
                            blocks = call_perm[[1]],
                            par_type = call_perm[[2]],
                            par_value = call_perm[[3]][b,],
                            perm = call_perm[[4]][b],
                            type = call_perm[[5]],
                            quiet = call_perm[[6]],
                            superblock = call_perm[[7]],
                            scheme = call_perm[[8]],
                            tol = call_perm[[9]], 
                            scale = call_perm[[10]], 
                            scale_block = call_perm[[11]], 
                            connection = call_perm[[12]], 
                            method = call_perm[[13]], 
                            bias = call_perm[[15]], 
                            init = call_perm[[16]], 
                            ncomp = call_perm[[17]], 
                            tau = call_perm[[18]], 
                            sparsity = call_perm[[18]]),
                          cl = cl)
    parallel::stopCluster(cl)
    rm("call_perm", envir = .GlobalEnv)
    }else{
      W = pbapply::pbsapply(seq(length(perm_parallel)), 
                            function(b) rgcca_permutation_k(
                              blocks = blocks,
                              par_type = par[[1]],
                              par_value = par_value_parallel[b,],
                              perm = perm_parallel[b],
                              type = type,
                              quiet = quiet,
                              superblock = superblock,
                              scheme = scheme,
                              tol = tol, 
                              scale = scale, 
                              scale_block = scale_block, 
                              connection = connection, 
                              method = method, 
                              bias = bias, 
                              init = init, 
                              ncomp = ncomp, 
                              tau = tau, 
                              sparsity = sparsity)
                            )
    }
    
    crits = W[!perm_parallel]
    permcrit = matrix(W[perm_parallel], nrow = nrow(par[[2]]), 
                      ncol = n_perms , byrow = TRUE)
    means = apply(permcrit, 1, mean, na.rm = T)
    sds = apply(permcrit, 1, sd, na.rm = T)
    par <- par[[2]]
    pvals <- sapply(seq(NROW(par)), function(k) mean(permcrit[k, ] >= crits[k]))
    zs <- sapply(seq(NROW(par)), function(k){ 
      z <- (crits[k] - mean(permcrit[k, ])) / (sd(permcrit[k, ]))
      if (is.na(z) || z == "Inf")
        z <- 0
      return(z)
    }
    )
    
    rownames(par) = 1:NROW(par)
    if(superblock){
      coln = c(names(blocks), "superblock")
    }else{
      coln = names(blocks)
    }
    colnames(par) = coln
    
    structure(list(call = call, zstat = zs, 
                   bestpenalties = par[which.max(zs), ],  
                   permcrit = permcrit, means = means, sds = sds, 
                   crit = crits, pvals = pvals, penalties = par, 
                   blocks = blocks), 
              class = "permutation"
              )
}
