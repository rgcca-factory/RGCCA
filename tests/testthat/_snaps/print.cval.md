# print_cval

    Code
      res <- rgcca_cv(blocks, response = 3, method = "rgcca", par_type = "tau",
        par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE)
    Output
      
    Code
      print(res, type = "quantile")
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1), response=3 
      There are J = 3 blocks.
      The design matrix is:
                  agriculture industry politic
      agriculture           0        0       1
      industry              0        0       1
      politic               1        1       0
      
      The factorial scheme is used.
      
      Tuning parameters (tau) used: 
        agriculture industry politic
      1           0      0.2     0.3
      2           0      0.0     0.0
      
      Validation: kfold with 5 folds and 1 run(s)) 
      Prediction model: lm 
      
        Tuning parameters Median error    Q1    Q3
      1       0.0/0.2/0.3        1.067 1.021 1.129
      2       0.0/0.0/0.0        0.984 0.947 1.274
      
      The best combination is: 0.0 0.2 0.3 for a mean CV error of 1.120 

---

    Code
      res <- rgcca_cv(blocks_classif, response = 3, method = "rgcca", par_type = "tau",
        par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE,
        prediction_model = "lda")
    Output
      
    Code
      print(res, type = "sd")
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1), response=3 
      There are J = 3 blocks.
      The design matrix is:
                  agriculture industry politic
      agriculture           0        0       1
      industry              0        0       1
      politic               1        1       0
      
      The factorial scheme is used.
      
      Tuning parameters (tau) used: 
        agriculture industry politic
      1           0      0.2       0
      2           0      0.0       0
      
      Validation: kfold with 5 folds and 1 run(s)) 
      Prediction model: lda 
      
        Tuning parameters Mean error Mean - Sd Mean + Sd
      1       0.0/0.2/0.0    0.08222  -0.00225   0.16670
      2       0.0/0.0/0.0    0.12889   0.01088   0.24690
      
      The best combination is: 0.0 0.2 0.0 for a mean CV error of  0.08222 

---

    Code
      res <- rgcca_cv(blocks, validation = "loo", response = 3, method = "sgcca",
        par_type = "sparsity", n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE)
    Output
      
    Code
      print(res, type = "sd")
    Output
      Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1), response=3 
      There are J = 3 blocks.
      The design matrix is:
                  agriculture industry politic
      agriculture           0        0       1
      industry              0        0       1
      politic               1        1       0
      
      The factorial scheme is used.
      
      Tuning parameters (sparsity) used: 
        agriculture industry politic
      1       1.000    1.000   1.000
      2       0.577    0.707   0.577
      
      Validation: loo 
      Prediction model: lm 
      
        Tuning parameters Mean error Mean - Sd Mean + Sd
      1    1.00/1.00/1.00      0.880     0.406     1.354
      2    0.58/0.71/0.58      0.932     0.438     1.426
      
      The best combination is: 1 1 1 for a mean CV error of 0.880 

