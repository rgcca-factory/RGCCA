data("Russett")
##############################################
# Test on the risk of having null variance   #
# variables in at least one bootstrap sample #
##############################################
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry    = Russett[, 4:5],
  politic     = Russett[, 6:11]
)

ncomp <- 1
# Rent is trapped.
blocks$agriculture$rent <- 0
blocks$agriculture$rent[1:4] <- 1
rgcca_out <- rgcca(blocks, ncomp = ncomp)

# When `pval = 1`, `generate_resampling` fails to identify `rent` as a risky
# variable, both when bootstraps are balanced or not.
set.seed(8882)
sample_out_balanced <- generate_resampling(
  rgcca_res = rgcca_out, n_boot = 4,
  balanced = T, pval = 1
)
sample_out_unbalanced <- generate_resampling(
  rgcca_res = rgcca_out, n_boot = 4,
  balanced = F, pval = 1
)
test_that("pVAL_high_noRiskyVAR", {
  expect_null(sample_out_balanced$sd_null)
  expect_null(sample_out_unbalanced$sd_null)
})

# Now, if `pval` is set to default, when `balanced = T` and
# `keep_all_variables = F`, a warning is generated to inform that variable `rent`
# is removed and `rent` is indeed removed.
set.seed(8882)
test_that("generate_resampling_missing_val_identification", {
  sample_out <- expect_warning(
    generate_resampling(rgcca_res = rgcca_out, n_boot = 4, balanced = T),
    paste0(
      "Variables:  rent appear to be of null ",
      "variance in some bootstrap samples and thus ",
      "were removed from all samples. \n",
      " ==> RGCCA is run again without these variables."
    )
  )
  expect_equal(
    names(sample_out$sd_null$agriculture),
    "rent"
  )
})

# Same situation, but this time, `pval` is set to its default value and it is
# specifically ask that all variables are kept. It is thus checked that `rent`
# is still there.
set.seed(8882)
sample_out <- generate_resampling(
  rgcca_res = rgcca_out, n_boot = 4,
  balanced = T, keep_all_variables = T
)
test_that("generate_resampling_keepAllVAriables", {
  expect_null(sample_out$sd_null)
})


#############################################
#   Test with 2 null variances variables    #
#############################################
# Now `rent` and `death` are trapped to be of null variance.
# Four tests are performed :
#   - "generate_resampling_NUL_variance_1" : when `balanced = T`
#      and `keep_all_variables = F`, a warning to inform that `rent` and `death`
#      are removed is raised.
#   - "generate_resampling_NUL_variance_2" : when `balanced = F`
#      and `keep_all_variables = F`, a warning to inform that `rent` and `death`
#      are removed is raised.
#   - "generate_resampling_NUL_variance_3" : when `balanced = T`
#      and `keep_all_variables = T`, an error is raised as it is impossible
#      to keep all variables here because some have null variances.
#   - "generate_resampling_NUL_variance_4" : when `balanced = T or F`
#      and `keep_all_variables = F`, check that `death` and `rent` are indeed
#      removed.
N <- NROW(Russett)
rgcca_out$call$raw$agriculture[, "rent"] <- rep(0, N)
rgcca_out$call$raw$politic[, "death"] <- rep(2, N)

test_that("generate_resampling_NUL_variance_1", {
  sample_out_balanced_1 <- expect_warning(
    generate_resampling(rgcca_res = rgcca_out, n_boot = 4, balanced = T),
    paste0(
      "Variables:  rent - death appear to be of null ",
      "variance in some bootstrap samples and thus ",
      "were removed from all samples. \n",
      " ==> RGCCA is run again without these variables."
    )
  )

  sample_out_balanced_2 <- expect_warning(
    generate_resampling(rgcca_res = rgcca_out, n_boot = 4, balanced = F),
    paste0(
      "Variables:  rent - death appear to be of null ",
      "variance in some bootstrap samples and thus ",
      "were removed from all samples. \n",
      " ==> RGCCA is run again without these variables."
    )
  )

  expect_equal(unname(unlist(sapply(
    sample_out_balanced_1$sd_null,
    function(x) names(x)
  ))), c("rent", "death"))
  expect_equal(unname(unlist(sapply(
    sample_out_balanced_2$sd_null,
    function(x) names(x)
  ))), c("rent", "death"))
})

test_that(
  "generate_resampling_NUL_variance_2",
  expect_error(
    generate_resampling(
      rgcca_res = rgcca_out, n_boot = 4,
      balanced = T, keep_all_variables = T
    ),
    paste0(
      "Impossible to define all bootstrap samples ",
      "without variables with null variance. Please ",
      "consider removing these variables:  rent - death.",
      " Please, consider unbalanced bootstrap by ",
      "setting 'balanced' to FALSE."
    )
  )
)

#########################################
#   Test with 2 very risky variables    #
#########################################
# Now `rent` and `death` are trapped to be very risky variables (only 1 observation
# differs from the others).
# Four tests are performed :
#   - "generate_resampling_veryRisky_1" : when `balanced = T`
#      and `keep_all_variables = F`, a warning to inform that `rent` and `death`
#      are removed is raised.
#   - "generate_resampling_veryRisky_2" : when `balanced = F`
#      and `keep_all_variables = F`, a warning to inform that `rent` and `death`
#      are removed is raised.
#   - "generate_resampling_veryRisky_3" : when `balanced = T`
#      and `keep_all_variables = T`, an error is raised as it is highly unlikely
#      to keep all variables here because some have almost null variances.
#   - "generate_resampling_veryRisky_4" : when `balanced = T or F`
#      and `keep_all_variables = F`, check that `death` and `rent` are indeed
#      removed.
#   - "generate_resampling_veryRisky_5" : when `balanced = F`
#      and `keep_all_variables = T`, check that no variable is removed.
N <- NROW(Russett)
rgcca_out$call$raw$agriculture[, "rent"] <- c(1, rep(0, N - 1))
rgcca_out$call$raw$politic[, "death"] <- c(1, rep(2, N - 1))
set.seed(553)
test_that("generate_resampling_veryRisky_1", {
  sample_out_balanced_1 <- expect_warning(
    generate_resampling(rgcca_res = rgcca_out, n_boot = 4, balanced = T),
    paste0(
      "Variables:  rent - death appear to be of null ",
      "variance in some bootstrap samples and thus ",
      "were removed from all samples. \n",
      " ==> RGCCA is run again without these variables."
    )
  )

  sample_out_balanced_2 <- expect_warning(
    generate_resampling(rgcca_res = rgcca_out, n_boot = 4, balanced = F),
    paste0(
      "Variables:  rent - death appear to be of null ",
      "variance in some bootstrap samples and thus ",
      "were removed from all samples. \n",
      " ==> RGCCA is run again without these variables."
    )
  )

  expect_equal(unname(unlist(sapply(
    sample_out_balanced_1$sd_null,
    function(x) names(x)
  ))), c("rent", "death"))
  expect_equal(unname(unlist(sapply(
    sample_out_balanced_2$sd_null,
    function(x) names(x)
  ))), c("rent", "death"))
})

set.seed(5553)
test_that(
  "generate_resampling_veryRisky_2",
  expect_error(
    generate_resampling(
      rgcca_res = rgcca_out, n_boot = 4,
      balanced = T, keep_all_variables = T
    ),
    paste0(
      "Impossible to define all bootstrap samples ",
      "without variables with null variance. Please ",
      "consider removing these variables:  rent - death.",
      " Please, consider unbalanced bootstrap by ",
      "setting 'balanced' to FALSE."
    )
  )
)

set.seed(53)
sample_out_balanced <- generate_resampling(
  rgcca_res = rgcca_out, n_boot = 4,
  balanced = F, keep_all_variables = T
)
test_that("generate_resampling_veryRisky_5", {
  expect_null(sample_out_balanced$sd_null)
})


##################################################################
#   Test with 2 very risky variables on a block of 2 variables   #
##################################################################
# Now `rent` and `death` are trapped to be very risky variables (only 1 observation
# differs from the others).
# Four tests are performed :
#   - "generate_resampling_ALL_Block_1" : when `balanced = T`
#      and `keep_all_variables = F`, an error is raised as it want to remove
#      all the variables from block `industry`.
#   - "generate_resampling_ALL_Block_2" : when `balanced = F`
#      and `keep_all_variables = F`, an error is raised as it want to remove
#      all the variables from block `industry`.
#   - "generate_resampling_ALL_Block_3" : when `balanced = T or F`
#      and `keep_all_variables = T` with a different random initialization,
#      no error is raised as no variable needs to be removed.
rgcca_out$call$raw$industry[, "gnpr"] <- c(1, rep(0, N - 1))
rgcca_out$call$raw$industry[, "labo"] <- c(1, rep(2, N - 1))
set.seed(54)
test_that(
  "generate_resampling_ALL_Block_1",
  expect_error(
    generate_resampling(
      rgcca_res = rgcca_out, n_boot = 4,
      balanced = T
    ),
    paste0(
      "The variance of all the variables from blocks:",
      "  industry appear to be null in some bootstrap ",
      "samples. Please consider removing them."
    )
  )
)
set.seed(52)
test_that(
  "generate_resampling_ALL_Block_2",
  expect_error(
    generate_resampling(
      rgcca_res = rgcca_out, n_boot = 4,
      balanced = F
    ),
    paste0(
      "The variance of all the variables from blocks:",
      "  industry appear to be null in some bootstrap ",
      "samples. Please consider removing them."
    )
  )
)

set.seed(1047)
sample_out_balanced <- generate_resampling(
  rgcca_res = rgcca_out, n_boot = 4,
  balanced = T
)
set.seed(6576)
sample_out_unbalanced <- generate_resampling(
  rgcca_res = rgcca_out, n_boot = 4,
  balanced = F
)
test_that("generate_resampling_ALL_Block_3", {
  expect_null(sample_out_balanced$sd_null)
  expect_null(sample_out_unbalanced$sd_null)
})
