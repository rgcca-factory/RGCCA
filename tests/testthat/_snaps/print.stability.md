# print.stability produces the expected text

    Code
      res <- rgcca_stability(fit.sgcca, n_boot = 10, verbose = FALSE)
      print(res)
    Output
      Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1),
      response=NULL, comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
                  agriculture industry politic
      agriculture           0        1       1
      industry              1        0       1
      politic               1        1       0
      
      The factorial scheme is used.
      Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 1.1576 
      
      The sparsity parameter used for agriculture is: 1 (with 2 variables selected)
      The sparsity parameter used for industry is: 1 (with 2 variables selected)
      The sparsity parameter used for politic is: 1 (with 3 variables selected)

