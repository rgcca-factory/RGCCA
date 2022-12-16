# print.bootstrap prints the expected string

    Code
      res <- rgcca_bootstrap(fit.rgcca, n_boot = 5, n_cores = 1, verbose = FALSE)
      print(res)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block='inertia', init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,2), response=NULL 
      There are J = 2 blocks.
      The design matrix is:
                  agriculture industry
      agriculture           0        1
      industry              1        0
      
      The factorial scheme is used.
      
      Extracted statistics on the block-weight vectors from 5 bootstrap samples 
      Component: 1 
             estimate       mean         sd lower_bound upper_bound bootstrap_ratio
      gini  0.6283955  0.6094169 0.07770200   0.5041928  0.69705679       8.0872498
      farm  0.7635057  0.7631563 0.06131671   0.7116233  0.85564322      12.4518377
      rent -0.1489231 -0.1229535 0.17046994  -0.3168363  0.06241138      -0.8736034
      gnpr  0.7397265  0.7489023 0.06103677   0.6888806  0.81202239      12.1193585
      labo -0.6729076 -0.6574286 0.07029069  -0.7248742 -0.58362603      -9.5732109
                pval adjust.pval
      gini 0.0000000   0.0000000
      farm 0.0000000   0.0000000
      rent 0.6666667   0.6666667
      gnpr 0.0000000   0.0000000
      labo 0.0000000   0.0000000
      Component: 2 
            estimate      mean         sd lower_bound upper_bound bootstrap_ratio
      gnpr 0.6729076 0.6574286 0.07029069   0.5836260   0.7248742        9.573211
      labo 0.7397265 0.7489023 0.06103677   0.6888806   0.8120224       12.119358
           pval adjust.pval
      gnpr    0           0
      labo    0           0

# print.bootstrap prints the expected string 2

    Code
      res <- rgcca_bootstrap(fit.rgcca, n_boot = 2, n_cores = 1, verbose = FALSE)
      print(res, type = "loadings")
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block='inertia', init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,2), response=NULL 
      There are J = 2 blocks.
      The design matrix is:
                  agriculture industry
      agriculture           0        1
      industry              1        0
      
      The factorial scheme is used.
      
      Extracted statistics on the block-loading vectors from 2 bootstrap samples 
      Component: 1 
             estimate       mean         sd lower_bound upper_bound bootstrap_ratio
      gini  0.9802706  0.9751404 0.11875857   0.9712584   0.9790224       19.404415
      farm  0.9784655  0.9741430 0.68204716   0.9559171   0.9923688        3.313861
      rent  0.3395018  0.3476065 0.12979257   0.2712455   0.4239675        2.723802
      gnpr  0.9571358  0.9586003 0.02730771   0.9571136   0.9600870       69.965580
      labo -0.9479563 -0.9463188 0.03029737  -0.9484440  -0.9441936      -59.781615
           pval adjust.pval
      gini    0           0
      farm    0           0
      rent    0           0
      gnpr    0           0
      labo    0           0
      Component: 2 
            estimate      mean          sd lower_bound upper_bound bootstrap_ratio
      gnpr 0.2896397 0.1129396 0.010890096   0.1057176   0.1201616        27.38019
      labo 0.3184005 0.4834443 0.007703263   0.4794791   0.4874095        42.82161
           pval adjust.pval
      gnpr    0           0
      labo    0           0

