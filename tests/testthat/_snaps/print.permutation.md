# print.permutation prints the expected text

    Code
      res <- rgcca_permutation(blocks, par_type = "tau", par_length = 2, n_perms = 5,
        n_cores = 1, verbose = FALSE)
      print(res)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1), response=NULL 
      There are J = 3 blocks.
      The design matrix is:
                  agriculture industry politic
      agriculture           0        1       1
      industry              1        0       1
      politic               1        1       0
      
      The factorial scheme is used.
      
      Tuning parameters (tau) used: 
        agriculture industry politic
      1           1        1       1
      2           0        0       0
      
        Tuning parameters Criterion Permuted criterion    sd zstat p-value
      1             1/1/1     0.717              0.150 0.058  9.77       0
      2             0/0/0     2.422              0.783 0.162 10.12       0
      
      The best combination is: 0, 0, 0 for a z score of 10.1 and a p-value of 0.

# print.permutation prints the expected text 2

    Code
      blocks2 <- rep(blocks, 3)
      names(blocks2) <- NULL
      res <- rgcca_permutation(blocks2, par_type = "ncomp", par_length = 2, n_perms = 2,
        n_cores = 1, verbose = FALSE)
      print(res)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1,1,1,1,1,1,1), response=NULL 
      There are J = 9 blocks.
      The design matrix is:
             block1 block2 block3 block4 block5 block6 block7 block8 block9
      block1      0      1      1      1      1      1      1      1      1
      block2      1      0      1      1      1      1      1      1      1
      block3      1      1      0      1      1      1      1      1      1
      block4      1      1      1      0      1      1      1      1      1
      block5      1      1      1      1      0      1      1      1      1
      block6      1      1      1      1      1      0      1      1      1
      block7      1      1      1      1      1      1      0      1      1
      block8      1      1      1      1      1      1      1      0      1
      block9      1      1      1      1      1      1      1      1      0
      
      The factorial scheme is used.
      
      Tuning parameters (ncomp) used: 
        block1 block2 block3 block4 block5 block6 block7 block8 block9
      1      2      2      2      2      2      2      2      2      2
      2      1      1      1      1      1      1      1      1      1
      
        Tuning parameters Criterion Permuted criterion    sd zstat p-value
      1             Set 1      16.6              1.032 0.144 107.6       0
      2             Set 2      15.6              0.953 0.258  56.6       0
      
      The best combination is: 2, 2, 2, 2, 2, 2, 2, 2, 2 for a z score of 108 and a p-value of 0.

