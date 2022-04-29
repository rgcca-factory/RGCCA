# print_cval

    Code
      res <- rgcca_cv(blocks, response = 3, method = "rgcca", par_type = "tau",
        par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE)
    Output
      
    Code
      print(res)
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
      1       0.0/0.2/0.3        0.772 0.709 0.874
      2       0.0/0.0/0.0        0.722 0.661 0.872
      
      The best combination is: 0 0 0 for a mean CV error of 0.722 

---

    Code
      res <- rgcca_cv(blocks, validation = "loo", response = 3, method = "sgcca",
        par_type = "sparsity", n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE)
    Output
      
    Code
      print(res, bars = "sd")
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
      2       0.577    0.707   0.408
      
      Validation: loo 
      Prediction model: lm 
      
        Tuning parameters Mean error Mean - Sd Mean + Sd
      1    1.00/1.00/1.00      0.606     0.359     0.853
      2    0.58/0.71/0.41      0.621     0.377     0.866
      
      The best combination is: 1 1 1 for a mean CV error of 0.606 

---

    Code
      blocks2 <- c(blocks, blocks)
      names(blocks2) <- NULL
      res <- rgcca_cv(blocks2, validation = "kfold", k = 2, response = 3, method = "rgcca",
        par_type = "tau", n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE)
    Output
      
    Code
      print(res, bars = "stderr")
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
      1 Tuning parameter set       0.692            0.682            0.702
      2 Tuning parameter set       0.808            0.687            0.929
      
      The best combination is: 1 1 1 1 1 1 for a mean CV error of 0.692 

