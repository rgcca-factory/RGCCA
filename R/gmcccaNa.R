gmgccaNa=function(blocks, method, connection = 1 - diag(length(blocks)),
                 tau = rep(1, length(blocks)),
                 ncomp = rep(1, length(blocks)),
                 scheme = "centroid", scale = TRUE,
                 init = "svd", bias = TRUE, tol = 1e-08,
                 verbose = TRUE, scale_block = TRUE, prescaling = FALSE,
                 quiet = FALSE, regularisation_matrices = NULL,
                 ranks = rep(1, length(blocks)), penalty_coef = 0)
{
	indNA=lapply(blocks, function(x){return(which(is.na(x), arr.ind = TRUE))})

  if(method=="complete") { A=intersection_list(blocks) }
  else if (is.function(method)) { A=method(blocks) }
	else { stop_rgcca("Only \"complete\" method is implemented to handle missing
	                  data for gMGCCA") }

  fit = gmgcca(A, C = connection, tau = tau, ncomp = ncomp,
              verbose = verbose, scale = scale, init = init, bias = bias,
              scale_block = scale_block, scheme = scheme,
              tol = tol, prescaling = prescaling, quiet = quiet,
              regularisation_matrices = regularisation_matrices,
              ranks = ranks, penalty_coef = penalty_coef)

 return(list(imputed_blocks = A, rgcca = fit, method, indNA = indNA))

}
