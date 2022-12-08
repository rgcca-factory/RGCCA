# print_cval

    Code
      res <- rgcca_cv(blocks, response = 3, method = "rgcca", par_type = "tau",
        par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE)
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
      
        Tuning parameters Median RMSE   Q1   Q3
      1       0.0/0.2/0.3        1.07 1.02 1.13
      2       0.0/0.0/0.0        1.09 1.03 1.13
      
      The best combination is: 0.0 0.2 0.3 for a mean RMSE of 1.12 

---

    Code
      res <- rgcca_cv(blocks_classif, response = 3, method = "rgcca", par_type = "tau",
        par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE,
        prediction_model = "lda")
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
      
        Tuning parameters Mean Accuracy Mean - Sd Mean + Sd
      1       0.0/0.2/0.0         0.893     0.787     0.999
      2       0.0/0.0/0.0         0.871     0.753     0.989
      
      The best combination is: 0.0 0.2 0.0 for a mean Accuracy of 0.893 

---

    Code
      res <- rgcca_cv(blocks, validation = "loo", metric = "MAE", response = 3,
        method = "sgcca", par_type = "sparsity", n_run = 1, n_cores = 1, par_length = 2,
        verbose = FALSE)
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
      1        1.00     1.00    1.00
      2        0.58     0.71    0.58
      
      Validation: loo 
      Prediction model: lm 
      
        Tuning parameters Mean MAE Mean - Sd Mean + Sd
      1    1.00/1.00/1.00    0.880     0.406     1.354
      2    0.58/0.71/0.58    0.930     0.437     1.422
      
      The best combination is: 1 1 1 for a mean MAE of 0.880 

---

    Code
      res <- rgcca_cv(blocks_classif, response = 3, method = "sgcca", par_type = "sparsity",
        n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE, prediction_model = "lda",
        metric = "Kappa")
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
      1        1.00     1.00       0
      2        0.58     0.71       0
      
      Validation: kfold with 5 folds and 1 run(s)) 
      Prediction model: lda 
      
        Tuning parameters Mean Kappa Mean - Sd Mean + Sd
      1    1.00/1.00/0.00      0.693     0.341     1.045
      2    0.58/0.71/0.00      0.618     0.172     1.064
      
      The best combination is: 1 1 0 for a mean Kappa of 0.693 

