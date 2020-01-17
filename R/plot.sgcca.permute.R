plot.sgcca.permute <- function (x) {
  c1s<- x$penalties
  ccs <- x$covs
  nperms <- x$nperms
  zstats <- x$zstat
  ccperms <- x$covperms
  par(mfrow = c(2, 1))
  plot(1:NROW(c1s), ccs, main = "Covariances For Real/Permuted Data", 
       xlab = "Index of Tuning Parameter Set", ylab = "Covariances", 
       ylim = range(ccperms, ccs))
  points(1:NROW(c1s), ccs, type = "l")
  for (i in 1:nperms) {
    points(1:NROW(c1s), ccperms[, i], col = "green")
  }
  plot(1:NROW(c1s), zstats, main = "Z", xlab = "Index of Tuning Parameter Set", 
       ylab = "Z score")
  lines(1:NROW(c1s), zstats)
}