#' Impute with RGCCA method
#'
#' This method is used for the implementation of EM algorithm for missing data
#'
#' @inheritParams select_analysis
#' @inheritParams rgccad
#' @param verbose FALSE by default. If TRUE, displays results of convergence for
#'  each iteration
#' @param threshold threshold to be reached to assess convergence.
#' @param reg 'y' by default. Reconstruction is made by regression of the blocks
#'  rows on w ('w'),the blocks columns on Y ('y')  or without any regression
#'  with y ('no')
#' @param ni An integer for the maximal number of iterations before convergence
#' @param naxis number of component to select for the estimation of missing data
#' @param scale TRUE if the variables should be standardized, FALSE ifelse
#' @param scale_block scaling of the blocks: "lambda1" for the first eigen value
#'  of the block (as in MFA)  or "inertia" for the sum of eigenvalues (when
#'  scale=TRUE, equivalent to divide by sqrt(p_j)).
#' @return \item{blocks}{blocks list of imputed matrices giving the \eqn{J}
#'  blocks of variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}}
#' @return \item{stab}{Convergence criterion : abs(1-obj_k/obj_{k-1})}
#' @return \item{obj}{Vector containing the mean square error between the
#'  predict values and the original non missing values at each iteration}
#' @return \item{crit}{RGCCA criterion}
#' @return \item{indNA}{Position of missing values}
#' @export
impute <-
  function(blocks,
           connection,
           tau,
           ni = 50,
           tol = 1e-08,
           ncomp = NULL,
           naxis = 1,
           scale = TRUE,
           scale_block = TRUE,
           scheme = "centroid",
           bias = TRUE,
           superblock = FALSE,
           verbose = FALSE,
           threshold = 0.001,
           reg = "y",
           quiet = FALSE) {
    # Initializations
    diff <- objective <- old <- criterion <- list()
    continue <- TRUE
    critRGCCA <- c()
    old[[1]] <- -Inf
    i <- 1
    nvar <- sapply(blocks, NCOL)
    nsuj <- dim(blocks[[1]])[1]
    J <- length(blocks)
    # Getting missing values pattern in the superblock
    indNA <- which(is.na(do.call(cbind, blocks)), arr.ind = TRUE)
    indNA2 <- lapply(blocks, function(x) {
      return(which(is.na(x), arr.ind = TRUE))
    })
    # Replacing NA by colmeans
    Alist <- impute_colmeans(blocks)

    while (continue) {
      diff[[i]] <- objective[[i]] <- old[[i + 1]] <- criterion[[i]] <- list()
      AlistScaled <- scaling(
        Alist,
        scale = scale,
        bias = bias,
        scale_block = scale_block
      )
      concatenedBlocks <- Reduce(cbind, Alist)
      fit.rgcca <- rgccad(
        blocks = AlistScaled,
        tau = tau,
        connection = connection,
        ncomp = ncomp,
        scheme = scheme,
        init = "svd",
        verbose = verbose,
        tol = tol
      )
      critRGCCA <- c(
        critRGCCA,
        fit.rgcca$crit[[1]][length(fit.rgcca$crit[[1]])]
      )
      # Reconstruction
      Xhat <- list()
      centeredXhat <- list()
      stdev <- moy <- sigma <- list()
      for (j in 1:J) {
        w <- fit.rgcca$astar[[j]]
        w <- w[, 1:naxis]
        # Getting back average and std
        moy[[j]] <- matrix(attr(AlistScaled[[j]], "scaled:center"),
          nrow = NROW(AlistScaled[[j]]),
          ncol = NCOL(AlistScaled[[j]]), byrow = TRUE
        )
        if (scale) {
          stdev[[j]] <- matrix(
            attr(AlistScaled[[j]], "scaled:scale"),
            nrow = NROW(AlistScaled[[j]]),
            ncol = NCOL(AlistScaled[[j]]),
            byrow = TRUE
          )
        } else {
          if (scale_block) {
            stdev[[j]] <- sqrt(NCOL(AlistScaled[[j]]))
          }
          if (!scale_block) {
            stdev[[j]] <- 1
          }
        }
        y <- fit.rgcca$Y[[j]][, naxis]
        centeredXhat[[j]] <- xhat_estimation(AlistScaled[[j]], y, naxis)
        residuals <- apply(AlistScaled[[j]], 2, function(x) {
          (lm(x ~ 0 + y)$residuals)
        })
        sigma[[j]] <- sqrt(sum(residuals^2 / (J * nsuj)))
        Xhat[[j]] <- (centeredXhat[[j]]) * stdev[[j]] + moy[[j]]
      }
      concatenedXhat <- Reduce(cbind, Xhat)
      concatenedBlocks[indNA] <- concatenedXhat[indNA]

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
      if (i > ni) {
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
        indNA = indNA2
      )
    )
  }
