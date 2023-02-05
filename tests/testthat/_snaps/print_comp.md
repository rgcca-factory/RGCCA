# print_comp prints the expected text

    Code
      fit.rgcca <- rgcca(blocks)
      res <- print_comp(fit.rgcca, n = 1, i = 1, outer = FALSE)
      cat(res)
    Output
      Comp. 1 (69.9%)

# print_comp prints the expected text 2

    Code
      fit.rgcca <- rgcca(blocks, ncomp = 2, sparsity = c(1, 1, 0.5))
      res <- print_comp(fit.rgcca, n = 2, i = 3, outer = FALSE)
      cat(res)
    Output
      Comp. 2 (3 variables, 18.2%)

# print_comp prints the expected text 3

    Code
      fit.rgcca <- rgcca(blocks)
      res <- print_comp(fit.rgcca, outer = TRUE)
      cat(res)
    Output
      First outer AVE: 60.2%

# print_comp prints the expected text 4

    Code
      fit.rgcca <- rgcca(blocks, ncomp = 4, superblock = TRUE)
      res <- print_comp(fit.rgcca, outer = TRUE)
      cat(res)
    Output
      First corrected outer AVE:  53.5% & 15.1%

