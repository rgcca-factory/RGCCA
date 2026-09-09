#' Bootstrap confidence intervals and p-values
#'
#' Bootstrap confidence intervals and p-values for evaluating the
#' significance/stability of the block-weight vectors produced by S/RGCCA.
#' @param rgcca_res A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}}).
#' @param n_boot The number of bootstrap samples (default: 100).
#' @param n_cores The number of cores used for parallelization.
#' @param verbose A logical value indicating if the progress of the bootstrap
#' procedure is reported.
#' @return A rgcca_bootstrap object that can be printed and plotted.
#' @return \item{n_boot}{The prinnumber of bootstrap samples, returned
#' for further use.}
#' @return \item{rgcca}{The RGCCA object fitted on the original data.}
#' @return \item{bootstrap}{A data.frame with the block weight vectors and
#' loadings computed on each bootstrap sample.}
#' @return \item{stats}{A data.frame of statistics summarizing the bootstrap
#' data.frame.}
#' @examples
#' # Bootstrap confidence intervals and p-values for RGCCA
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:8]
#' )
#'
#' fit_rgcca <- rgcca(blocks, ncomp = 1)
#'
#' boot_out <- rgcca_bootstrap(fit_rgcca, n_boot = 20, n_cores = 1,
#'                             verbose = TRUE)
#'
#' print(boot_out)
#' plot(boot_out, type = "weight", block = 1:3, comp = 1,
#'      display_order = FALSE)
#'
#'
#' \dontrun{
#'
#'  # Download the dataset's package at http://biodev.cea.fr/sgcca/ and install
#'  # it from the package archive file.
#'  # You can do it with the following R commands:
#'  if (!("gliomaData" %in% rownames(installed.packages()))) {
#'    destfile <- tempfile()
#'    download.file(
#'      "http://biodev.cea.fr/sgcca/gliomaData_0.4.tar.gz", destfile
#'    )
#'    install.packages(destfile, repos = NULL, type = "source")
#'  }
#'
#'  data("ge_cgh_locIGR", package = "gliomaData")
#'  blocks <- ge_cgh_locIGR$multiblocks
#'  Loc <- factor(ge_cgh_locIGR$y)
#'  levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#'  blocks [[3]] <- Loc
#'
#'
#'  fit_sgcca <- rgcca(blocks, response = 3,
#'     sparsity = c(.071, .2, 1), ncomp = 1,
#'    scheme = "factorial",
#'    verbose = TRUE
#'  )
#'
#'  print(fit_sgcca)
#'
#'  boot_out <- rgcca_bootstrap(fit_sgcca, n_boot = 50, n_cores = 2)

#'  plot(boot_out, block = 1:2, type = "weight",
#'        comp = 1, n_mark = 300000,
#'        display_order = FALSE)
#' }
#' @export
#' @seealso \code{\link[RGCCA]{plot.rgcca_bootstrap}},
#' \code{\link[RGCCA]{summary.rgcca_bootstrap}}
rgcca_bootstrap <- function(rgcca_res, n_boot = 100,
                            n_cores = 1, verbose = TRUE) {

 stability <- is(rgcca_res, "rgcca_stability")
  if (stability) {
    message(
      "All the parameters were imported from the fitted rgcca_stability",
      " object."
    )
    rgcca_res <- rgcca_res$rgcca_res
  }
    extract <- function(A, .dim, .value) {
    idx.list <- lapply(dim(A), seq_len)
    idx.list[[.dim]] <- .value
    do.call(`[`, c(list(A), idx.list))
}

  # If sparse model, we perform bootstrap only on the selected variables
  if (!stability && tolower(rgcca_res$call$method) %in% sparse_methods()) {
    if (verbose) {
      message(
        "Only selected variables were used for bootstrapping. see ",
        "rgcca_stability()."
      )
    }

    # Remove superblock variables from keep_var as the superblock is generated
    # from the kept variables
    J <- length(rgcca_res$call$blocks)
   
    for (i in 1:(length(rgcca_res$call$blocks)-1)){
    
      if (length(dim(rgcca_res$call$blocks[[i]]))>2){
        for (m in 1:(length(dim(rgcca_res$call$blocks[[i]]))-1)){
          if (!is.null(dim(rgcca_res$factors[[i]][m][[1]]))){
            aux <- rgcca_res$factors[[i]][m][[1]][,1]
          }
          else{
            aux <- rgcca_res$factors[[i]][m][[1]]
          }
         
         names(aux)=seq_len(length(aux))
     
        x=unlist(aux)!=0
        
        if (m==1){
      keep=c(unlist(aux)[x])

        }
        else{
          keep = c(keep,unlist(aux)[x])
          
        }
       # print(rgcca_res$a)
      
        
          rgcca_res$call$blocks[[i]]=extract(rgcca_res$call$blocks[[i]],m+1,seq_len(length(aux))[x])
         


    } 
 }
     else{
      aux <- rgcca_res$a[[i]]
      if (!is.null(dim(aux))){
        aux <- aux[, 1]
      }
      keepvar=unlist(aux)!=0
     
       rgcca_res$call$blocks[[i]]=rgcca_res$call$blocks[[i]][,keepvar,drop=FALSE]

     }}
  #  keep_var <- lapply(
  #    rgcca_res$a[-(J + 1)],
  #    function(x) unique(which(x != 0, arr.ind = TRUE)[, 1])
  #  )
  #  if (rgcca_res$opt$disjunction) {
  #    keep_var[[rgcca_res$call$response]] <- 1
  #  }

    #rgcca_res$call$blocks <- Map(
    #  function(x, y) x[, y, drop = FALSE], rgcca_res$call$blocks, keep_var
   # )
    rgcca_res$call$tau <-
      rgcca_res$call$sparsity <- rep(1, length(rgcca_res$blocks))
    
    rgcca_res <- rgcca(rgcca_res,rgcca_res$call$method)
  }
  

  check_integer("n_boot", n_boot)


  ### Create bootstrap samples
  # If there is a disjunctive response block, sample bootstrap samples
  # with a stratified strategy
  if (rgcca_res$opt$disjunction) {
    folds <- caret::createFolds(
      rep(rgcca_res$call$blocks[[rgcca_res$call$response]][, 1], n_boot),
      k = n_boot, list = TRUE,
      returnTrain = FALSE
    )
    idx <- rep(seq_len(NROW(rgcca_res$call$blocks[[1]])), n_boot)
    v_inds <- lapply(folds, function(f) idx[f])
  } else {
    v_inds <- lapply(seq_len(n_boot), function(i) {
      sample(seq_len(NROW(rgcca_res$call$blocks[[1]])), replace = TRUE)
    })
  }
  
  ### Run RGCCA on the bootstrap samples
  W <- par_pblapply(v_inds, function(b) {
    rgcca_bootstrap_k(
      rgcca_res = rgcca_res,
      inds = b
    )
  }, n_cores = n_cores, verbose = verbose)
  
  factors_df <- NULL
  
  multi_blocks <- sapply(rgcca_res$call$blocks, function(x) is.array(x) && length(dim(x)) > 2)
  
  if (any(multi_blocks)) {
    factors_df <- do.call(rbind, lapply(seq_along(W), function(b) { 
      do.call(rbind, lapply(seq_along(W[[b]]$F), function(j) {
        F_list <- W[[b]]$F[[j]] 
        if (is.null(F_list)) return(NULL) 
        do.call(rbind, lapply(seq_along(F_list), function(k) {
          M <- F_list[[k]] 
          n <- nrow(M)
          p <- ncol(M)
          if (is.null(p)){
            p=1
            vars <- names(M)}

          else{
    
          vars <- rownames(M)}
          
          data.frame(
            boot  = b,
            block = names(rgcca_res$a)[j],
            mode  = k,
            var   = rep(vars, times = p),
            comp  = rep(seq_len(p), each = n),
            value = as.vector(M),
            row.names = NULL
          )
        }))
      }))
    }))
    if (nrow(factors_df) == 0) factors_df <- NULL
  }
 
  # Test unimodality for each mode using the dip test
  # Since we are in a multivariate setting, we apply the dip test on the
  # first principal component.
  W2 <- lapply(W, function(x) { x$F <- NULL; x })
  res <- format_bootstrap_list(W2, rgcca_res)
  res <- check_sign_comp(rgcca_res, res)


  res_f <- check_sign_comp_factors(rgcca_res, factors_df)
  
  # common statistics
  if (!is.null(res_f) && nrow(res_f) > 0) {
    res_f$type <- "factors"
    cols_communes <- intersect(colnames(res), colnames(res_f))
    res_sub <- res[, cols_communes, drop = FALSE]
    resf_sub <- res_f[, cols_communes, drop = FALSE]
    res_glob <- rbind(res_sub, resf_sub)
  } else {
    res_glob <- res
  }
  
  ### Extract statistics from the results of the bootstrap
  stats <- rgcca_bootstrap_stats(res_glob, rgcca_res, length(W))

  return(structure(list(n_boot     = n_boot,
                        rgcca      = rgcca_res,
                        bootstrap  = res, 
                        factors    = res_f,
                        stats      = data.frame(stats)),
                   class = "rgcca_bootstrap"))
}