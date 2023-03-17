gmgcca_aux_varNa=function(blocks, method, connection = 1 - diag(length(blocks)),
                  tau = rep(1, length(blocks)),
                  ncomp = rep(1, length(blocks)),
                  scheme = "centroid", scale = TRUE,
                  init = "svd", bias = TRUE, tol = 1e-08,
                  verbose = TRUE, scale_block = TRUE, prescaling = FALSE,
                  quiet = FALSE, ranks = rep(1, length(blocks)),
                  n_run = 1, n_cores = 1)
{
  indNA=lapply(blocks, function(x){return(which(is.na(x), arr.ind = TRUE))})

  if(method=="complete") { A=intersection_list(blocks) }
  else if (is.function(method)) { A=method(blocks) }
  else { stop_rgcca("Only \"complete\" method is implemented to handle missing
	                  data for gMGCCA") }

  fit = gmgcca_aux_var(A, C = connection, tau = tau, ncomp = ncomp,
               verbose = verbose, scale = scale, init = init, bias = bias,
               scale_block = scale_block, scheme = scheme,
               tol = tol, prescaling = prescaling, quiet = quiet,
               ranks = ranks)

  return(list(imputed_blocks = A, rgcca = fit, method, indNA = indNA))

}
