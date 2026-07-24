#' Center and scale a list of blocks
#' @inheritParams rgcca
#' @param na.rm A logical, if TRUE, NA values are replaced by 0 to
#' compute scaling parameters.
#' @noRd

scaling <- function(blocks, scale = TRUE, bias = TRUE,
                    scale_block = "inertia", na.rm = TRUE) {
                       
  if (isTRUE(scale_block)) scale_block <- "inertia"
  sqrt_N <- sqrt(NROW(blocks[[1]]) + bias - 1)

  blocks <- lapply(blocks, function(x) {
    # Store dim and dimnames
    dim_x <- dim(x)
    dimnames_x <- dimnames(x)
    if (is.numeric(scale)){
      
      if (length(dim(x))>2){
          dim(x)=matrix(c(dim(x)[1], dim(x)[scale+1], prod(dim(x)[-c(1,scale+1)])),nrow=1)

      for (i in 1:dim(x)[2]){
     
            x[1:dim(x)[1], i, 1:dim(x)[3]] <- scale2(x[1:dim(x)[1],i,1:dim(x)[3]], scale=TRUE,  bias=bias)
      }
      scale=TRUE
         
            
      x <- matrix(x, nrow = nrow(x))}
      else{
      x<- scale2(x, scale = TRUE, bias = bias)


      }}
    else{
      ##baseline
      if (length(x)>2){

        x <- matrix(x, nrow = nrow(x))
      }

        x<- scale2(x, scale = TRUE, bias = bias)
    }
          
       # Scale each block by a constant if requested
    if (scale_block == "lambda1") {
      x <- scale_lambda1(x, sqrt_N, scale, na.rm = na.rm)
    } else if (scale_block == "inertia") {
      x <- scale_inertia(x, sqrt_N, scale, na.rm = na.rm)
    }
      dim(x)=dim_x
   

    # Go back to a tensor
    y <- array(x, dim = dim_x)
    dimnames(y) <- dimnames_x
    attr(y, "scaled:center") <- attr(x, "scaled:center")
    attr(y, "scaled:scale") <- attr(x, "scaled:scale")

    return(y)
  })

  return(blocks)
}
  