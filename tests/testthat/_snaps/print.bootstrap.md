# print.bootstrap prints the expected string

    Code
      res <- rgcca_bootstrap(fit.rgcca, n_boot = 5, n_cores = 1, verbose = FALSE)
      print(res)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,2),
      response=NULL, comp_orth=TRUE 
      There are J = 2 blocks.
      The design matrix is:
                  agriculture industry
      agriculture           0        1
      industry              1        0
      
      The factorial scheme is used.
      
      Extracted statistics from 5 bootstrap samples.
      Block-weight vectors for component 1: 
           estimate   mean     sd lower_bound upper_bound bootstrap_ratio  pval
      gini    0.628  0.609 0.0777       0.504      0.6971           8.087 0.000
      farm    0.764  0.763 0.0613       0.712      0.8556          12.452 0.000
      rent   -0.149 -0.123 0.1705      -0.317      0.0624          -0.874 0.667
      gnpr    0.740  0.749 0.0610       0.689      0.8120          12.119 0.000
      labo   -0.673 -0.657 0.0703      -0.725     -0.5836          -9.573 0.000
           adjust.pval
      gini       0.000
      farm       0.000
      rent       0.667
      gnpr       0.000
      labo       0.000

# print.bootstrap prints the expected string 2

    Code
      res <- rgcca_bootstrap(fit.rgcca, n_boot = 2, n_cores = 1, verbose = FALSE)
      print(res, type = "loadings", comp = 2, block = 2)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block='inertia',
      init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,2),
      response=NULL, comp_orth=TRUE 
      There are J = 2 blocks.
      The design matrix is:
                  agriculture industry
      agriculture           0        1
      industry              1        0
      
      The factorial scheme is used.
      
      Extracted statistics from 2 bootstrap samples.
      Block-loading vectors for component 2: 
           estimate  mean     sd lower_bound upper_bound bootstrap_ratio pval
      gnpr    0.290 0.113 0.0109       0.106       0.120            27.4    0
      labo    0.318 0.483 0.0077       0.479       0.487            42.8    0
           adjust.pval
      gnpr           0
      labo           0

