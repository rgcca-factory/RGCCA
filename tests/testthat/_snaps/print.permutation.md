# print.permutation

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
      
        Tuning parameters Criterion Permuted criterion sd     zstat  p-value
      1 1/1/1              0.717     0.150              0.058  9.770  0.000 
      2 0/0/0              2.422     0.783              0.162 10.121  0.000 
      
      The best combination is: 0, 0, 0 for a z score of 10.1 and a p-value of 0.

---

    Code
      blocks2 <- c(blocks, blocks)
      names(blocks2) <- NULL
      res <- rgcca_permutation(blocks2, par_type = "ncomp", par_length = 2, n_perms = 2,
        n_cores = 1, verbose = FALSE)
      print(res)
    Output
      Call: method='rgcca', superblock=FALSE, scale=TRUE, scale_block=TRUE, init='svd', bias=TRUE, tol=1e-08, NA_method='nipals', ncomp=c(1,1,1,1,1,1), response=NULL 
      There are J = 6 blocks.
      The design matrix is:
             block1 block2 block3 block4 block5 block6
      block1      0      1      1      1      1      1
      block2      1      0      1      1      1      1
      block3      1      1      0      1      1      1
      block4      1      1      1      0      1      1
      block5      1      1      1      1      0      1
      block6      1      1      1      1      1      0
      
      The factorial scheme is used.
      
      Tuning parameters (ncomp) used: 
        block1 block2 block3 block4 block5 block6
      1      2      2      2      2      2      2
      2      1      1      1      1      1      1
      
        Tuning parameters      Criterion Permuted criterion sd       zstat   
      1 Tuning parameter set 1   6.2353    0.3610             0.0365 161.0061
      2 Tuning parameter set 2   5.8807    0.5198             0.0297 180.6396
        p-value 
      1   0.0000
      2   0.0000
      
      The best combination is: 1, 1, 1, 1, 1, 1 for a z score of 181 and a p-value of 0.

