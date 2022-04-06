#' Imputes with superblock method
#' @inheritParams impute
#' @export
impute_sb <- function(blocks,
                      tau = rep(1, length(blocks) + 1),
                      ni = 5,
                      tol = 1e-08,
                      ncomp = NULL,
                      naxis = 1,
                      scale = TRUE,
                      scale_block = TRUE,
                      scheme = "centroid",
                      bias = TRUE,
                      verbose = FALSE,
                      threshold = 1e-6,
                      reg = "y",
                      quiet = FALSE) {
  # Getting missing values pattern in the superblock
  indNA <- which(is.na(do.call(cbind, blocks)), arr.ind = TRUE)
  indNA2 <- lapply(blocks, function(x) {
    return(which(is.na(x), arr.ind = TRUE))
  })
  # Initializations
  nvar <- sapply(blocks, NCOL)
  nsuj <- dim(blocks[[1]])[1]
  J <- length(blocks)
  diff <- objective <- old <- criterion <- list()
  continue <- TRUE
  critRGCCA <- c()
  old[[1]] <- -Inf
  i <- 1
  connection <- matrix(0, J + 1, J + 1)
  connection[1:J, J + 1] <- 1
  connection[J + 1, 1:J] <- 1
  ncomp <- c(rep(1, J), naxis)

  # Imputing blocks by the colmeans
  Alist <- impute_colmeans(blocks)
  
  
  while (continue) {
    diff[[i]] <- objective[[i]] <- old[[i + 1]] <- criterion[[i]] <- list()
    # scaling original blocks
    AlistScaled0 <- scaling(
      Alist,
      scale = scale,
      bias = bias,
      scale_block = scale_block
    )
    centerVec <- scaleVec <- c()

    for (j in 1:J) {
      centerVec <- c(centerVec, attr(AlistScaled0[[j]], "scaled:center"))
      scaleVec <- c(scaleVec, attr(AlistScaled0[[j]], "scaled:scale"))
    }

    concatenedBlocks <- Reduce(cbind, Alist)
    # superblock
    scaledConcatenedBlocks <- Reduce(cbind, AlistScaled0)

    # list including the superblock
    AlistScaled <- AlistScaled0
    AlistScaled[[J + 1]] <- scaledConcatenedBlocks
    attr(AlistScaled[[J + 1]], "scaled:center") <- centerVec
    attr(AlistScaled[[J + 1]], "scaled:scale") <- scaleVec
    # run rgcca on the obtained list
    fit.rgcca <- rgccad(
      blocks = AlistScaled,
      tau = tau,
      connection = connection,
      ncomp = ncomp,
      scheme = scheme,
      init = "svd",
      verbose = verbose,
      tol = tol,
      superblock = TRUE
    )
    critRGCCA <- c(critRGCCA, fit.rgcca$crit[[1]][length(fit.rgcca$crit[[1]])])
    # Reconstruction
    Xhat <- list()
    stdev <- moy <- sigma <- list()
    moy[[J + 1]] <-
      matrix(
        attr(AlistScaled[[J + 1]], "scaled:center"),
        nrow = NROW(AlistScaled[[J + 1]]),
        ncol = NCOL(AlistScaled[[J + 1]]),
        byrow = TRUE
      )
    stdev[[J + 1]] <- matrix(
      attr(AlistScaled[[J + 1]], "scaled:scale"),
      nrow = NROW(AlistScaled[[J + 1]]),
      ncol = NCOL(AlistScaled[[J + 1]]),
      byrow = TRUE
    )
    
    y <- fit.rgcca$Y[[J + 1]][, 1:naxis]
    centeredXhat <- xhat_estimation(scaled_superblock=AlistScaled[[J+1]], y=y, naxis=naxis)
    residuals <- apply(AlistScaled[[J + 1]], 2, function(x) {
      (lm(x ~ 0 + y)$residuals)
    })
    sigma[[J + 1]] <- sqrt(sum(residuals^2 / (J * nsuj)))
    print(stdev[[J+1]][1,])
    print(moy[[J+1]][1,])
    Xhat[[J + 1]] <- (centeredXhat) * stdev[[J + 1]]+ moy[[J + 1]]
    concatenedXhat <- Reduce(cbind, Xhat)
    concatenedBlocks[indNA] <- concatenedXhat[indNA]
    scaledConcatenedBlocks <- scale2(concatenedBlocks,
      scale = scale, bias = FALSE
    )

    # Criteria of convergence are calculated
    #----------------------------------------
    diff[[i]] <- concatenedBlocks - concatenedXhat
    diff[[i]][indNA] <- 0
    objective[[i]] <-
      sqrt(sum(diff[[i]]^2) / (
        dim(concatenedBlocks)[1] *
          dim(concatenedBlocks)[2] - sum(is.na(concatenedBlocks))
      ))
    criterion[[i]] <- abs(1 - objective[[i]] / old[[i]])
    old[[i + 1]] <- objective[[i]]
    if (!is.nan(criterion[[i]])) {
      if (criterion[[i]] < threshold && (i > 2)) {
        continue <- FALSE
      }
      if (objective[[i]] < threshold) {
        continue <- FALSE
        if (verbose) {
          cat(
            "The algorithm converged to a stationary point after",
            i - 1,
            "iterations \n"
          )
        }
      }
    }
    if (i >= ni) {
      continue <- FALSE
      if (!quiet) {
        warning(
          paste(
            "The RGCCA imputation algorithm did not converge after ",
            ni,
            " iterations"
          )
        )
      }
    }


    # Reaffecting Alist
    #----------------------
    Alist <- superblock_to_matrix(concatenedBlocks, nvar)
    names(Alist) <- names(blocks)
    i <- i + 1
  }
  return(
    list(
      A = Alist,
      crit = unlist(critRGCCA),
      stab = unlist(criterion),
      obj = unlist(objective),
      stdev = stdev,
      sigma = sigma,
      moy = moy,
      indNA = indNA2
    )
  )
}
