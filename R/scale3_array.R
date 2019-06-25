scale3_array<-function (A, center = TRUE, scale = TRUE, bias = TRUE, method_scale = "3-way") 
{
  if (center == TRUE & scale == TRUE) {
    if (method_scale == "3-way"){
      B   = t(apply(A, 1, c))
      B   = scale(B, center = TRUE, scale = FALSE)
      C   = reshape3(B, n_col = dim(A)[2])
      std = sqrt(apply(C, 2, function(x) cov2(x, bias = bias)))
      
      if (any(std==0)) {
        sprintf("there were %d constant variables",sum(std==0))
        std[std==0]=1
      }
      
      B = B/matrix(rep(std, NROW(B)), NROW(B), NCOL(B), byrow = TRUE)
      A = array(c(as.matrix(B)), dim(A), dimnames = dimnames(A))
      
      names(std)               = dimnames(A)[[2]]
      attr(A, "scaled:scale")  = std
      center                   = matrix(attr(B, "scaled:center"), nrow = dim(A)[2], ncol = dim(A)[3])
      rownames(center)         = dimnames(A)[[2]]
      colnames(center)         = dimnames(A)[[3]]
      attr(A, "scaled:center") = center
      return(A)
    }else{
      B   = t(apply(A, 1, c))
      B   = scale(B, center = TRUE, scale = FALSE)
      std = sqrt(apply(B, 2, function(x) cov2(x, bias = bias)))
      
      if (any(std==0)) {
        sprintf("there were %d constant variables",sum(std==0))
        std[std==0]=1
      }
      
      B = B/matrix(rep(std, NROW(B)), NROW(B), NCOL(B), byrow = TRUE)
      A = array(c(as.matrix(B)), dim(A), dimnames = dimnames(A))
      
      std                      = matrix(std, nrow = dim(A)[2], ncol = dim(A)[3])
      rownames(std)            = dimnames(A)[[2]]
      colnames(std)            = dimnames(A)[[3]]
      attr(A, "scaled:scale")  = std
      center                   = matrix(attr(B, "scaled:center"), nrow = dim(A)[2], ncol = dim(A)[3])
      rownames(center)         = dimnames(A)[[2]]
      colnames(center)         = dimnames(A)[[3]]
      attr(A, "scaled:center") = center
      return(A)
    }
   }
  if (center == TRUE & scale == FALSE) {
    B   = t(apply(A, 1, c))
    B   = scale(B, center = TRUE, scale = FALSE)
    A   = array(c(as.matrix(B)), dim(A), dimnames = dimnames(A))
    
    center                   = matrix(attr(B, "scaled:center"), nrow = dim(A)[2], ncol = dim(A)[3])
    rownames(center)         = dimnames(A)[[2]]
    colnames(center)         = dimnames(A)[[3]]
    attr(A, "scaled:center") = center
    return(A)
  }
  if (center == FALSE & scale == TRUE) {
    if (method_scale == "3-way"){
      B   = t(apply(A, 1, c))
      C   = reshape3(B, n_col = dim(A)[2])
      std = sqrt(apply(C, 2, function(x) cov2(x, bias = bias)))
      
      if (any(std==0)) {
        sprintf("there were %d constant variables",sum(std==0))
        std[std==0]=1
      }
      
      B = B/matrix(rep(std, NROW(B)), NROW(B), NCOL(B), byrow = TRUE)
      A = array(c(as.matrix(B)), dim(A), dimnames = dimnames(A))
      
      names(std)              = dimnames(A)[[2]]
      attr(A, "scaled:scale") = std
      return(A)
    }else{
      B   = t(apply(A, 1, c))
      std = sqrt(apply(B, 2, function(x) cov2(x, bias = bias)))
      
      if (any(std==0)) {
        sprintf("there were %d constant variables",sum(std==0))
        std[std==0]=1
      }
      
      B = B/matrix(rep(std, NROW(B)), NROW(B), NCOL(B), byrow = TRUE)
      A = array(c(as.matrix(B)), dim(A), dimnames = dimnames(A))
      
      std                      = matrix(std, nrow = dim(A)[2], ncol = dim(A)[3])
      rownames(std)            = dimnames(A)[[2]]
      colnames(std)            = dimnames(A)[[3]]
      attr(A, "scaled:scale")  = std
      return(A)
    }
  }
}
