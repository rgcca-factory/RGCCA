# print_call prints the expected text

    Code
      res <- rgcca(blocks, ncomp = 1, tau = 1, scheme = "horst", connection = 1 -
        diag(3), scale = FALSE, scale_block = "lambda1", superblock = FALSE,
      response = NULL, NA_method = "na.omit")
      print_call(res$call)
    Output
      Call: method='rgcca', superblock=FALSE, scale=FALSE, scale_block='lambda1',
      init='svd', bias=TRUE, tol=1e-08, NA_method='na.omit', ncomp=c(1,1,1),
      response=NULL, comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
             block1 block2 block3
      block1      0      1      1
      block2      1      0      1
      block3      1      1      0
      
      The horst scheme is used.

# print_call prints the expected text 2

    Code
      res <- rgcca(blocks, ncomp = 2, tau = 1, scheme = "centroid", scale = TRUE,
        scale_block = TRUE, superblock = FALSE, response = 3, NA_method = "na.ignore")
      print_call(res$call)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd',
      bias=TRUE, tol=1e-08, NA_method='na.ignore', ncomp=c(2,2,2), response=3,
      comp_orth=TRUE 
      There are J = 3 blocks.
      The design matrix is:
             block1 block2 block3
      block1      0      0      1
      block2      0      0      1
      block3      1      1      0
      
      The centroid scheme is used.

