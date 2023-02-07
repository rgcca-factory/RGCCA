# print.rgcca produces the expected text

    Code
      res <- rgcca(blocks, ncomp = 1, tau = c(1, 0.7, 0.9))
      print(res)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1),
      response=NULL, comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
             block1 block2 block3
      block1      0      1      1
      block2      1      0      1
      block3      1      1      0
      
      The factorial scheme is used.
      Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 1.1932 
      
      The regularization parameter used for block1 is: 1
      The regularization parameter used for block2 is: 0.7
      The regularization parameter used for block3 is: 0.9

# print.rgcca produces the expected text 2

    Code
      tau <- matrix(c(1, 1, 1, 0, 0, 0), nrow = 2, byrow = TRUE)
      res <- rgcca(blocks, ncomp = 2, tau = tau)
      print(res)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(2,2,2),
      response=NULL, comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
             block1 block2 block3
      block1      0      1      1
      block2      1      0      1
      block3      1      1      0
      
      The factorial scheme is used.
      Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 1.7093 
      
      The regularization parameters used are: 
           block1 block2 block3
      [1,]      1      1      1
      [2,]      0      0      0

# print.rgcca produces the expected text 3

    Code
      res <- rgcca(blocks, ncomp = 1, sparsity = c(0.7, 1, 0.9), method = "sgcca")
      print(res)
    Output
      Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1),
      response=NULL, comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
             block1 block2 block3
      block1      0      1      1
      block2      1      0      1
      block3      1      1      0
      
      The factorial scheme is used.
      Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.9348 
      
      The sparsity parameter used for block1 is: 0.7 (with 2 variables selected)
      The sparsity parameter used for block2 is: 1 (with 2 variables selected)
      The sparsity parameter used for block3 is: 0.9 (with 2 variables selected)

# print.rgcca produces the expected text 4

    Code
      sparsity <- matrix(c(0.8, 1, 0.9, 1, 0.9, 0.8), nrow = 2, byrow = TRUE)
      res <- rgcca(blocks, ncomp = 2, sparsity = sparsity, method = "sgcca")
      print(res)
    Output
      Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(2,2,2),
      response=NULL, comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
             block1 block2 block3
      block1      0      1      1
      block2      1      0      1
      block3      1      1      0
      
      The factorial scheme is used.
      Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.9958 
      
      The sparsity parameters used are: 
           block1 block2 block3
      [1,]    0.8    1.0    0.9
      [2,]    1.0    0.9    0.8
      The number of selected variables are: 
           block1 block2 block3
      [1,]      2      2      2
      [2,]      3      2      2

# print.rgcca produces the expected text 5

    Code
      res <- rgcca(blocks2, ncomp = 1, sparsity = c(0.7, 1, 0.9), response = 3)
      print(res)
    Output
      Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1),
      response=3, comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
             block1 block2 block3
      block1      0      0      1
      block2      0      0      1
      block3      1      1      0
      
      The factorial scheme is used.
      Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 1.0814 
      
      The sparsity parameter used for block1 is: 0.7 (with 2 variables selected)
      The sparsity parameter used for block2 is: 1 (with 2 variables selected)
      The regularization parameter used for block3 is: 0

