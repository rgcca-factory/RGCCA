xhat_estimation <- function(scaled_superblock, y, naxis) {
  if (naxis == 1) {
    gamma <- apply(scaled_superblock, 2, function(x) {
      lm(x ~ 0 + y)$coefficients[1]
    })
    centeredXhat <- matrix(y, ncol = naxis) %*%
      matrix(gamma, nrow = naxis)
  }
  if (naxis > 1) {
    gamma <- NULL
    for (k in 1:naxis) {
      # regression of the relevant block on the y
      gamma <- cbind(
        gamma,
        apply(
          scaled_superblock,
          2,
          function(x) lm(x ~ y[, k])$coefficients[2]
        )
      )
    }
    centeredXhat <- matrix(y, ncol = naxis) %*% t(gamma)
  }
  return(centeredXhat)
}
