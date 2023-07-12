gmgcca_pddNa=function(blocks, method, connection = 1 - diag(length(blocks)),
                              tau = rep(1, length(blocks)),
                              ncomp = rep(1, length(blocks)),
                              scheme = "centroid", scale = TRUE,
                              init = "svd", bias = TRUE, tol = 1e-08,
                              verbose = TRUE, scale_block = TRUE, prescaling = FALSE,
                              quiet = FALSE, n_run = 1, n_cores = 1, penalty_coef = 1)
{
  indNA=lapply(blocks, function(x){return(which(is.na(x), arr.ind = TRUE))})
  na.rm = FALSE
  if(method == "complete") A = intersection_list(blocks)
  if(method == "nipals"){na.rm = TRUE ; A = blocks}

  fit = gmgcca_pdd(A, connection, tau = tau, ncomp = ncomp,
                           verbose = verbose, scale = scale, init = init, bias = bias,
                           scale_block = scale_block, scheme = scheme,
                           tol = tol, prescaling = prescaling, quiet = quiet,
                           penalty_coef = penalty_coef)

  return(list(imputed_blocks = A, rgcca = fit, method, indNA = indNA))

}
