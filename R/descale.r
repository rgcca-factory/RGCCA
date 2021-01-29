descale=function(blocks){
res = lapply(blocks, 
             function(x){
               if(is.null(attr(x, "scaled:center"))) return(x) 
               else{
                 center_vec = attr(x, "scaled:center") 
                 scaling_vec = attr(x, "scaled:scale") 
                 center_mat = matrix(rep(center_vec, nrow(x)),
                                     ncol = length(center_vec), byrow = T)
                 scaling_mat = matrix(rep(scaling_vec, nrow(x)),
                                      ncol = length(center_vec), byrow = T) 
                 attr(x, "scaled:scale") <- NULL 
                 attr(x, "scaled:center") <- NULL 
                 y = x * scaling_mat + center_mat 
                 return(y)
               }
             }
            )
return(res)
}