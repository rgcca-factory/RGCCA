#' Available methods for RGCCA
#'
#' List the methods that can be used with the rgcca function.
#' @return A vector of the methods implemented with the rgcca function.
#' @examples
#' available_methods()
#' @export
available_methods <- function() {
  c(
    "rgcca", "sgcca", "pca", "spca", "pls", "spls", "cca",
    "ifa", "ra", "gcca", "maxvar", "maxvar-b", "maxvar-a",
    "mfa", "mcia", "mcoa", "cpca-1", "cpca-2", "cpca-4", "hpca", "maxbet-b",
    "maxbet", "maxdiff-b", "maxdiff", "sabscor",
    "ssqcor", "ssqcov-1", "ssqcov-2", "ssqcov",
    "sumcor", "sumcov-1", "sumcov-2", "sumcov", "sabscov-1",
    "sabscov-2", "netsgcca"
  )
}

one_block_methods <- function() {
  c("pca", "spca")
}

two_block_methods <- function() {
  c("cca", "ra", "ifa", "pls", "spls")
}

superblock_methods <- function() {
  c(
    "pca", "spca", "gcca", "maxvar", "maxvar-b", "maxvar-a",
    "mfa", "mcia", "mcoa", "cpca-1", "cpca-2", "cpca-4", "hpca"
  )
}

cov_methods <- function() {
  c(
    "pca", "spca", "pls", "ifa", "maxvar-a", "mfa", "mcia", "mcoa", "cpca-1",
    "cpca-2", "cpca-4", "hpca", "maxbet-b", "maxbet", "maxdiff-b",
    "maxdiff", "ssqcov-1", "ssqcov-2", "ssqcov", "sumcov-1",
    "sumcov-2", "sumcov", "sabscov-1", "sabscov-2"
  )
}

cor_methods <- function() {
  c("cca", "gcca", "maxvar", "maxvar-b", "sabscor", "ssqcor", "sumcor")
}

horst_methods <- function() {
  c(
    "pca", "spca", "pls", "spls", "cca", "ifa", "ra", "cpca-1",
    "maxbet", "maxdiff", "sumcor", "sumcov-1", "sumcov-2",
    "sumcov"
  )
}

factorial_methods <- function() {
  c(
    "gcca", "maxvar", "maxvar-b", "maxvar-a", "mfa", "mcia", "mcoa",
    "cpca-2", "maxbet-b", "maxdiff-b", "ssqcor", "ssqcor",
    "ssqcov-1", "ssqcov-2", "ssqcov"
  )
}

centroid_methods <- function() {
  c("sabscor", "sabscov-1", "sabscov-2")
}

x4_methods <- function() {
  c("cpca-4", "hpca")
}

sparse_methods <- function() {
  c("sgcca", "spls", "spca", "netsgcca")
}

generic_methods <- function() {
  c("rgcca", "sgcca", "netsgcca")
}
