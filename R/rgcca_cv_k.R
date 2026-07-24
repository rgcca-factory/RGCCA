#' Cross-validation
#'
#' Uses cross-validation to evaluate predictive model of RGCCA
#' @inheritParams rgcca_predict
#' @inheritParams rgcca
#' @inheritParams rgcca_bootstrap
#' @noRd
rgcca_cv_k <- function(rgcca_args, inds,prediction_model,params=NULL,tuning=NULL,
                       par_type, par_value, metric, upsample=FALSE,...) {
  
 
  rgcca_args[[par_type]] <- par_value
 #rgcca_args[['sparse_lambda']] <- par_value2
  blocks <- rgcca_args[["blocks"]]
  ind_pos_=setdiff(c(1:dim(rgcca_args[["blocks"]][[1]])[1]),inds)
 
 # intermediate <- lapply(
 #   blocks, function(x) subset_block_rows(x, ind_pos, drop = FALSE)
 # )
  if (upsample){
          ind_pos=NULL
          ind_pos$ind=ind_pos_
          
          ind_pos$response=rgcca_args[["blocks"]]
          ind_pos_new=upSample(y=as.factor(rgcca_args[["blocks"]]$response[ind_pos$ind]),
          x=ind_pos$ind)              
  }else{
    ind_pos_new=NULL
    ind_pos_new$x=ind_pos_
  }
  
  rgcca_args[["blocks"]] <- lapply(
    blocks, function(x) subset_block_rows(x, ind_pos_new$x, drop = FALSE)
  )

  # Fit RGCCA on the training blocks
  for (i in 1:length(rgcca_args[["blocks"]])){
      row.names(rgcca_args[["blocks"]][[i]]) <- NULL

  }
 
  
  res <- do.call(rgcca, rgcca_args)

  # Evaluate RGCCA on the validation blocks
  blocks_test <- lapply(
    blocks, function(x) subset_block_rows(x, inds, drop = FALSE)
  )
  names(blocks_test) <- names(res$blocks)

  return(rgcca_predict(
    res,
    metric = metric,
    blocks_test = blocks_test,
    prediction_model = prediction_model,params=params,tuning=tuning,
    ...
  )$score)
}
