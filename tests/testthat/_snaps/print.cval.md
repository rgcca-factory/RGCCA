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
      
        Tuning parameters Median error  2.5% 97.5%
      1       0.0/0.2/0.3        1.067 1.007 1.328
      2       0.0/0.0/0.0        1.143 0.893 1.392
      
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
      1       0.0/0.2/0.0     0.1289    0.0408    0.2169
      2       0.0/0.0/0.0     0.1067    0.0357    0.1776
      
      The best combination is: 0 0 0 for a mean CV error of 0.1067 

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

---

    Code
      blocks2 <- c(blocks, blocks)
      names(blocks2) <- NULL
      res <- rgcca_cv(blocks2, validation = "kfold", k = 2, response = 3, method = "rgcca",
        par_type = "tau", n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE)
    Output
      
    Code
      print(res, type = "stderr")
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1,1,1,1), response=3 
      There are J = 6 blocks.
      The design matrix is:
           [,1] [,2] [,3] [,4] [,5] [,6]
      [1,]    0    0    1    0    0    0
      [2,]    0    0    1    0    0    0
      [3,]    1    1    0    1    1    1
      [4,]    0    0    1    0    0    0
      [5,]    0    0    1    0    0    0
      [6,]    0    0    1    0    0    0
      
      The factorial scheme is used.
      
      Tuning parameters (tau) used: 
        [,1] [,2] [,3] [,4] [,5] [,6]
      1    1    1    1    1    1    1
      2    0    0    0    0    0    0
      
      Validation: kfold with 2 folds and 1 run(s)) 
      Prediction model: lm 
      
            Tuning parameters Mean error Mean - Std Error Mean + Std Error
      1 Tuning parameter set       0.775            0.673            0.878
      2 Tuning parameter set       0.657            0.494            0.821
      
      The best combination is: 0 0 0 0 0 0 for a mean CV error of 0.657 

