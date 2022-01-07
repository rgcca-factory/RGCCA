list_khatri_rao <- function(factors) {
  Reduce("khatri_rao", rev(factors))
}

kron_sum <- function(factors) {
  apply(list_khatri_rao(factors), 1, sum)
}

weighted_kron_sum <- function(factors, weights) {
  list_khatri_rao(factors) %*% weights
}

mix_weighted_kron_sum <- function(factors, weights, ncomp, ranks) {
  sapply(1:ncomp, function(k) {
    weighted_kron_sum(
      lapply(factors, "[", seq((k - 1) * ranks + 1, k * ranks)),
      weights[k, ]
    )
  })
}

kron_prod_q <- function(factors, mode, q) {
  D = length(factors)
  Reduce("%x%", rev(lapply(1:D, function(d) {
    if (mode == d) {
      return(diag(nrow(factors[[mode]])))
    } else {
      return(factors[[d]][, q])
    }
  })))
}

inv_sqrtm = function(M){
  eig        = eigen(M)
  M_inv_sqrt = eig$vectors %*% diag(eig$values^(-1/2), nrow = nrow(M)) %*% t(eig$vectors)
  return(M_inv_sqrt)
}

construct_projector = function(Mq, Mqmq, Wmq, tol = 1e-8) {
  # Start by finding the active orthogonality constraints
  tmp = Mqmq %*% Wmq
  idx = sapply(1:ncol(tmp), function(x) {
    norm(tmp[, x], type = "2") > tol
  })
  Wmq = Wmq[, idx, drop = FALSE]
  # Construct projector to satisfy orthogonality constraints
  if (dim(Wmq)[2] == 0) pi = 0
  else {
    tmp  = Mq %*% Mqmq %*% Wmq
    pi   = tmp %*% ginv(crossprod(tmp)) %*% t(tmp)
  }
  return(pi)
}

project_factor_q <- function(factors, mode, q, Mq, Mqmq) {
  Wmq = list_khatri_rao(lapply(factors, function(x) x[, 1:(q - 1), drop = F]))
  pi = construct_projector(Mq, Mqmq, Wmq)
  return(Mq %*% (diag(nrow(Mq)) - pi) %*% Mq %*% factors[[mode]][, q])
}

weighted_factor <- function(u, d, rank) {
  if (rank == 1)
    return(u)
  return((u %*% diag(d[1:rank])) / sqrt(sum(d[1:rank] ^ 2)))
}

alt_prod = function(XtX, factor, LEN, d, r, side = "left") {
  unfolded_dim = dim(XtX)
  folded_dim   = sapply(factor, nrow)
  if (side == "left") {
    XtX    = array(XtX, dim = c(folded_dim, unfolded_dim[2]))
    offset = 0
  }
  else {
    XtX    = array(XtX, dim = c(unfolded_dim[1], folded_dim))
    offset = 1
  }
  for (m in 1:(LEN - 1)) {
    if (d != m) {
      XtX = mode_product(XtX, factor[[m]][, r], mode = m + offset)
      folded_dim[m] = 1
    }
  }
  if (side == "left")
    return(drop(array(XtX, dim = c(folded_dim, unfolded_dim[2]))))
  else
    return(drop(array(XtX, dim = c(unfolded_dim[1], folded_dim))))
}
