#' Identify the most stable variables with SGCCA
#'
#' This function can be used to identify the most stable variables
#' identified as relevant by SGCCA. A Variable Importance in the Projection
#' (VIP) based criterion is used to identify the most stable variables.
#'
#' @inheritParams rgcca_bootstrap
#' @param keep A numeric vector indicating the proportion of variables per
#' block to select.
#' @param verbose A logical value indicating if the progress of the procedure
#' is reported.
#' @return A rgcca_stability object that can be printed and plotted.
#' @return \item{top}{A data.frame giving the indicator (VIP)
#' on which the variables are ranked.}
#' @return \item{n_boot}{The number of bootstrap samples, returned
#' for further use.}
#' @return \item{keepVar}{The indices of the most stable variables.}
#' @return \item{bootstrap}{A data.frame with the block weight vectors
#' computed on each bootstrap sample.}
#' @return \item{rgcca_res}{An RGCCA object fitted on the most stable
#' variables.}
#' @examples
#' \dontrun{
#'  ###########################
#'  # stability and bootstrap #
#'  ###########################
#'
#'  data("ge_cgh_locIGR", package = "gliomaData")
#'  blocks <- ge_cgh_locIGR$multiblocks
#'  Loc <- factor(ge_cgh_locIGR$y)
#'  levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#'  blocks[[3]] <- Loc
#'
#'  fit_sgcca <- rgcca(blocks,
#'     sparsity = c(.071, .2, 1),
#'     ncomp = c(1, 1, 1),
#'     scheme = "centroid",
#'     verbose = TRUE, response = 3
#' )
#'
#'  boot_out <- rgcca_bootstrap(fit_sgcca, n_boot = 100, n_cores = 1)
#'
#'  fit_stab <- rgcca_stability(fit_sgcca,
#'    keep = sapply(fit_sgcca$a, function(x) mean(x != 0)),
#'    n_cores = 1, n_boot = 10,
#'    verbose = TRUE
#'  )
#'
#'  boot_out <- rgcca_bootstrap(
#'    fit_stab, n_boot = 500, n_cores = 1, verbose = TRUE
#'  )
#'
#'  plot(boot_out, block = 1:2, n_mark = 2000, display_order = FALSE)
#' }
#' @export
rgcca_stability <- function(rgcca_res,
                            keep = vapply(
                              rgcca_res$a, function(x) mean(x != 0),
                              FUN.VALUE = 1.0
                            ),method=NULL,
                            n_boot = 100,
                            n_cores = 1,
                            verbose = TRUE) {
  stopifnot(tolower(rgcca_res$call$method) %in% sparse_methods())
  check_integer("n_boot", n_boot)
  check_integer("n_cores", n_cores, min = 0)
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
  #
  #if (method=="stgcca"){
  #  v_inds <- lapply(seq_len(n_boot), function(i) {
  #  sample(seq_len(NCOL(rgcca_res$call$blocks[[1]])), replace = TRUE)
  #})
  #}
    
  W <- par_pblapply(v_inds, function(b) {
    rgcca_bootstrap_k(
      rgcca_res = rgcca_res,
      inds = b, type = "AVE"
    )
  }, n_cores = n_cores, verbose = verbose)
  factors_df <- NULL
  if (method!='rgcca'){
  
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
          vars <- rownames(M)
          if (is.null(p)){
            p=1
          }
          if (is.null(n)){
            n=1
          }
          if (is.null(vars)){
            vars="selected"
          }
        

          
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
  }}
  
  # Test unimodality for each mode using the dip test
  # Since we are in a multivariate setting, we apply the dip test on the
  # first principal component.
  W2 <- lapply(W, function(x) { x$F <- NULL; x })
  
  res <- format_bootstrap_list(W2, rgcca_res)
  res <- check_sign_comp(rgcca_res, res)
  #res <- format_bootstrap_list(W, rgcca_res)
  if (method!="rgcca"){
  res_f<- check_sign_comp_factors(rgcca_res, factors_df)}
  
  if (method!="rgcca"){
  if (!is.null(res_f) && nrow(res_f) > 0) {
    res_f$type <- "factors"
    cols_communes <- intersect(colnames(res), colnames(res_f))
    res_sub <- res[, cols_communes, drop = FALSE]
    resf_sub <- res_f[, cols_communes, drop = FALSE]
    res_glob <- rbind(res_sub, resf_sub)
  } else {
    res_glob <- res
  }}else{
    res_glob=res
  }
    res_or=res

  res=res_glob
  J <- length(rgcca_res$blocks)

  if (rgcca_res$call$superblock == TRUE) {
    res <- res[res$block != names(rgcca_res$blocks)[J], ]
    rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[-J]
    rgcca_res$call$blocks <- rgcca_res$call$blocks[-J]
  }

  if (rgcca_res$opt$disjunction) {
     res <- res[res$block != names(rgcca_res$blocks)[rgcca_res$call$response], ]
     rgcca_res$AVE$AVE_X <- rgcca_res$AVE$AVE_X[-rgcca_res$call$response]
  }
  res_AVE <- res[res$type == "loadings", ]
  res <- res[res$type == "factors", ]

  # Compute var2block to later retrieve "block" from "var"
  var2block <- subset(res, res$comp == 1 & res$boot == 1)[, c("var", "block")]
  rownames(var2block) <- var2block$var
  var2block$var <- NULL
 

  
  grp <- c('col1')

# calculating mean of col2 based on col1 group
  #factors_df=factors_df %>% 
  #group_by(across(c("boot","comp","block"))) %>% 
  #summarize(value = mean(value))


  res$scores <-  res$value^2 *factors_df$value#
  
  top <- tapply(
    res$scores, list(var = res$var), mean
  )
  top <- cbind(top = top, block = var2block[names(top), ])
  top<- top[!is.na(top[,"block"]),]
  top_factors=top
  perc <- elongate_arg(keep, top)


  if (is.null(dim(rgcca_res$call$sparsity))) {
    if (rgcca_res$call$method!='stgcca'){
    if (rgcca_res$call$superblock == TRUE) {
      rgcca_res$call$sparsity <- rgcca_res$call$sparsity[-J]
    }
    perc[which(rgcca_res$call$sparsity == 1)] <- 1
  } }else {
    if (rgcca_res$call$method!='stgcca'){
    if (rgcca_res$call$superblock == TRUE) {
      rgcca_res$call$sparsity <- rgcca_res$call$sparsity[, -J]
    }
    perc[which(rgcca_res$call$sparsity[1, ] == 1)] <- 1
  }}

  # Keep a percentage of the variables with the top intensities
  keepVar<- lapply(seq_along(rgcca_res$AVE$AVE_X), function(j) {
    x <- top[top[, "block"] == j, "top"]
    order(x, decreasing = TRUE)[seq(round(perc[j] * length(x)))]
  })
  if (rgcca_res$opt$disjunction) {
    keepVar[[rgcca_res$call$response]] <- 1
  }
  res_AVE <- res_or[res_or$type != "weights", ]
  res <- res_or[res_or$type == "weights", ]


  # Compute var2block to later retrieve "block" from "var"
  var2block <- subset(res, res$comp == 1 & res$boot == 1)[, c("var", "block")]
  rownames(var2block) <- var2block$var
  var2block$var <- NULL
 



  res$scores <-  res$value^2 *res_AVE$value#
  
  top <- tapply(
    res$scores, list(var = res$var), mean
  )
  top <- cbind(top = top, block = var2block[names(top), ])
  top<- top[!is.na(top[,"block"]),]
  perc <- elongate_arg(keep, top)


  if (is.null(dim(rgcca_res$call$sparsity))) {
    if (rgcca_res$call$method!='stgcca'){
    if (rgcca_res$call$superblock == TRUE) {
      rgcca_res$call$sparsity <- rgcca_res$call$sparsity[-J]
    }
    perc[which(rgcca_res$call$sparsity == 1)] <- 1
  } }else {
    if (rgcca_res$call$method!='stgcca'){
    if (rgcca_res$call$superblock == TRUE) {
      rgcca_res$call$sparsity <- rgcca_res$call$sparsity[, -J]
    }
    perc[which(rgcca_res$call$sparsity[1, ] == 1)] <- 1
  }}
   extract <- function(A, .dim, .value) {
    idx.list <- lapply(dim(A), seq_len)
    idx.list[[.dim]] <- .value
    do.call(`[`, c(list(A), idx.list))
}
    

  # Keep a percentage of the variables with the top intensities
  keepVar2 <- lapply(seq_along(rgcca_res$AVE$AVE_X), function(j) {
    x <- top[top[, "block"] == j, "top"]
    order(x, decreasing = TRUE)[seq(round(perc[j] * length(x)))]
  })
  if (rgcca_res$opt$disjunction) {
    keepVar2[[rgcca_res$call$response]] <- 1
  }
 for (i in 1:(length(rgcca_res$call$blocks)-1)){
  if (length(dim(rgcca_res$call$blocks[[i]]))>2){
    for (m in 1:(length(dim(rgcca_res$call$blocks[[i]]))-1)){
     
      if (rgcca_res$call$sparsity[[i]][m]<1){
        perc= mean(rgcca_res$factors[[i]][[m]]!=0)
        
        x=top_factors[top_factors[, "block"] == i, "top"][unlist(dimnames(rgcca_res$call$blocks[[i]])[m+1])]
        if (m==1){
      keep=c(order(x, decreasing = TRUE)[seq(round(perc * length(x)))])

        }
        else{
          keep = list(keep,c(order(x, decreasing = TRUE)[seq(round(perc * length(x)))] ))
          
        }
      
            rgcca_res$call$blocks[[i]]=extract(rgcca_res$call$blocks[[i]],m+1,c(order(x, decreasing = TRUE)[seq(round(perc * length(x)))] ))


   
      } 

    }
  
    #rgcca_res$call$blocks[[i]]=rgcca_res$call$blocks[[i]][1:dim(rgcca_res$call$blocks[[i]])[1],keep[[1]],keep[[2]]]

  }
     else{
       rgcca_res$call$blocks[[i]]=rgcca_res$call$blocks[[i]][,keepVar2[[i]],drop=FALSE]

     }}
     
      


     



  original=rgcca_res$call$blocks
    



  rgcca_res$call$tau <-
    rgcca_res$call$sparsity <- rep(1, length(rgcca_res$call$blocks))


  rgcca_res <- rgcca(rgcca_res)




  return(structure(list(
    top = top,
    n_boot = n_boot,
    keepVar = keepVar,
    bootstrap = res,
    rgcca_res = rgcca_res
  ),
  class = "rgcca_stability"
  ))
}
