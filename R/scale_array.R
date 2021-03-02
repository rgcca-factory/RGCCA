# Standardization (to zero means and unit variances or mean squares) of arrays.
# In case of nD tensors (with n > 2), scaling is done "within" the first mode.
# As it may not be the most appropriate choice, the user should preprocess his
# data beforehand.
# @inheritParams rgccaNa
# @param A A numeric matrix.
# @param center A logical value. If center = TRUE, each column is translated to
# have zero mean.
# @return \item{A}{The centered and/or scaled matrix. The centering and
# scaling values (if any) are returned as attributes "scaled:center" and
# "scaled:scale".}
# @title Scaling and Centering of arrays.
# @export scale_array

scale_array <- function (A, center = TRUE, scale = TRUE, bias = TRUE)
{
  DIM = dim(A)
  if (length(DIM) > 2) {
    B   = matrix(as.vector(A), nrow = DIM[1]) # Mode 1 matricization of A

    if (scale == TRUE) {
      ms  = apply(B, 1, function(x) sqrt(sum(x ^ 2)))
      if (any(ms == 0)) {
        sprintf("there were %d constant variables", sum(ms == 0))
        ms[ms == 0] = 1
      }
      B   = apply(B, -1, function(x) x / ms)
    }

    B   = scale(B, center = center, scale = FALSE)
    A   = array(as.vector(B), dim(A), dimnames = dimnames(A))
    attr(A, "scaled:center") = attr(B, "scaled:center")
    if (scale == TRUE) {
      attr(A, "scaled:scale") = ms
    }
    return(A)
  } else {
    return(scale2(A, center = center, scale = scale, bias = bias))
  }
}

# scale_array <- function (A, center = TRUE, scale = TRUE, bias = TRUE)
# {
#   DIM = dim(A)
#   if (length(DIM) > 2) {
#     B   = matrix(as.vector(A), nrow = DIM[1]) # Mode 1 matricization of A
#     B   = scale(B, center = center, scale = FALSE)
#     C   = unfold(array(B, dim = DIM), mode = 2)
#
#     if (scale == TRUE) {
#       ms  = apply(C, 1, function(x) sqrt(sum(x ^ 2)))
#       if (any(ms == 0)) {
#         sprintf("there were %d constant variables", sum(ms == 0))
#         ms[ms == 0] = 1
#       }
#       C   = apply(C, -1, function(x) x / ms)
#     }
#
#     A   = array(as.vector(C), dim(A), dimnames = dimnames(A))
#     attr(A, "scaled:center") = attr(B, "scaled:center")
#     if (scale == TRUE) {
#       attr(A, "scaled:scale") = ms
#     }
#     return(A)
#   } else {
#     return(scale2(A, center = center, scale = scale, bias = bias))
#   }
# }
