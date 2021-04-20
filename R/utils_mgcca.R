list_khatri_rao <- function(factors) {
  Reduce("khatri_rao", rev(factors))
}

kron_sum <- function(factors) {
  apply(list_khatri_rao(factors), 1, sum)
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
  M_inv_sqrt = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
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
