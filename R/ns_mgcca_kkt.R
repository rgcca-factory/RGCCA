ns_mgcca_kkt = function(A, A_m, a, factors, weights, XtX, Y, g, dg, C,
                        ranks = rep(1, length(A)), bias = T) {

  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  n      <- DIM[[1]][1]
  Z      <- matrix(0, n, J)

  # Gradient of constraints on lambda
  dh_l_j = function(weights) {
    matrix(2 * weights)
  }

  # Gradient of norm constraints on w_jmr
  dk_l_jr = function(factors, XtX, rank, LEN, m, r) {
    Mqmq = alt_prod(XtX, factors, LEN, m, r, side = "left")
    Mq   = alt_prod(Mqmq, factors, LEN, m, r, side = "right")
    2 * Mq %*% factors[[m]][, r]
  }

  # Gradient of orthogonality constraints on w_jmr and w_jms
  dq_l_jrs = function(factors, XtX, rank, LEN, m, r, s) {
    ws   = Reduce("%x%", rev(lapply(factors, function(x) x[, s])))
    Mqmq = alt_prod(XtX, factors, LEN, m, r, side = "left")
    Mqmq %*% ws
  }

  # Gradient of the objective function
  grad_l = function(A_m, Y, dg, C, j) {
    dgx    = dg(cov2(Y[, j], Y, bias = bias))
    dgx    = matrix(rep(dgx, n), n, J, byrow = TRUE)
    Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
    res    = t(A_m[[j]]) %*% Z[, j]
  }

  df_w_jmr = function(grad, factors, weights, m, r) {
    weights[r] * t(kron_prod_q(factors, m, r)) %*% grad
  }

  df_l_j = function(grad, factors) {
    t(list_khatri_rao(factors)) %*% grad
  }

  sizes = c()
  for (j in 1:J) sizes = c(sizes, ranks[j])
  for (j in 1:J) sizes = c(sizes, ranks[j] * sum(DIM[[j]][-1]))
  total_size = sum(sizes)

  dh = list()
  df = matrix(0, nrow = total_size, ncol = 1)

  ### Constraints on lambda
  for (j in 1:J) {
    dh[[j]] = matrix(0, nrow = total_size, ncol = 1)
    idx     = (1 + sum(sizes[c(0:(j - 1))])):sum(sizes[c(0:j)])
    dh[[j]][idx, ] = dh_l_j(weights[[j]])
  }

  ### Norm constraints
  for (j in 1:J) {
    for (r in 1:ranks[j]) {
      dh[[length(dh) + 1]] = matrix(0, nrow = total_size, ncol = 1)
      idx  = (1 + sum(sizes[c(0:(J + j - 1))])):sum(sizes[c(0:(J + j))])
      dh_l = c()
      for (s in 1:ranks[j]) {
        if (s != r) {
          dh_l = rbind(dh_l, matrix(0, nrow = sum(DIM[[j]][-1]), ncol = 1))
        } else {
          for (m in 1:(LEN[j] - 1)) {
            dh_l = rbind(dh_l, dk_l_jr(factors[[j]], XtX[[j]], ranks[j], LEN[j], m, r))
          }
        }
      }
      dh[[length(dh)]][idx, ] = dh_l
    }
  }

  ### Orthogonality constraints
  for (j in 1:J) {
    for (r in 1:(ranks[j] - 1)) {
      for (s in (r + 1):ranks[j]) {
        dh[[length(dh) + 1]] = matrix(0, nrow = total_size, ncol = 1)
        idx  = (1 + sum(sizes[c(0:(J + j - 1))])):sum(sizes[c(0:(J + j))])
        dh_l = c()
        for (t in 1:ranks[j]) {
          if (t != r && t != s) {
            dh_l = rbind(dh_l, matrix(0, nrow = sum(DIM[[j]][-1]), ncol = 1))
          } else if (t == r) {
            for (m in 1:(LEN[j] - 1)) {
              dh_l = rbind(dh_l, dq_l_jrs(factors[[j]], XtX[[j]], ranks[j], LEN[j], m, r, s))
            }
          } else {
            for (m in 1:(LEN[j] - 1)) {
              dh_l = rbind(dh_l, dq_l_jrs(factors[[j]], XtX[[j]], ranks[j], LEN[j], m, s, r))
            }
          }
        }
        dh[[length(dh)]][idx, ] = dh_l
      }
    }
  }

  ### Gradient with respect to lambda
  for (j in 1:J) {
    idx  = (1 + sum(sizes[c(0:(j - 1))])):sum(sizes[c(0:j)])
    grad = grad_l(A_m, Y, dg, C, j)
    df[idx, ] = df_l_j(grad, factors[[j]])
  }

  ### Gradient with respect to w_jmr
  for (j in 1:J) {
    idx  = (1 + sum(sizes[c(0:(J + j - 1))])):sum(sizes[c(0:(J + j))])
    df_l = c()
    grad = grad_l(A_m, Y, dg, C, j)
    for (r in 1:ranks[j]) {
      for (m in 1:(LEN[j] - 1)) {
        df_l = rbind(df_l, df_w_jmr(grad, factors[[j]], weights[[j]], m, r))
      }
    }
    df[idx, ] = df_l
  }

  ### Put constraints together
  dh = do.call(cbind, dh)

  Pi = dh %*% ginv(crossprod(dh)) %*% t(dh)

  kkt = norm((diag(nrow(Pi)) - Pi) %*% df, type = "2") / norm(df, type = "2")
  return(kkt)
}
