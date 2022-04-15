#' The function sgccak() is called by sgcca() and does not have to be used by
#' the user. sgccak() enables the computation of SGCCA block components, outer
#' weight vectors, etc., for each block and each deflation stage.
#' @inheritParams select_analysis
#' @inheritParams rgccad
#' @inheritParams sgcca
#' @inheritParams rgccak
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{crit}{A vector of integer that contains for each component
#' the values of the analysis criteria across iterations.}
#' @return \item{converg}{Speed of convergence of the alogrithm to reach the
#' tolerance.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @title Internal function for computing the SGCCA parameters (SGCCA block
#' components, outer weight vectors etc.)
#' @importFrom Deriv Deriv
sgccak <- function(A, C, sparsity = rep(1, length(A)),
                   scheme = "centroid", tol = 1e-08,
                   init = "svd", bias = TRUE, verbose = FALSE,
                   quiet = FALSE, na.rm = TRUE, n_iter_max = 1000) {
  if (is.function(scheme)) {
    g <- scheme
    # check for parity of g
    ctrl <- !any(g(-5:5) != g(5:-5))
  } else {
    switch(scheme,
      "horst" = {
        g <- function(x) x
        ctrl <- FALSE
      },
      "factorial" = {
        g <- function(x) x^2
        ctrl <- TRUE
      },
      "centroid" = {
        g <- function(x) abs(x)
        ctrl <- TRUE
      }
    )
  }

  dg <- Deriv::Deriv(g, env = parent.frame())

  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- sapply(A, NCOL)
  const <- sparsity * sqrt(pjs)

  tmp <- sgcca_init(A, init, bias, na.rm, sparsity, pjs, J, n)
  a <- tmp$a
  Y <- tmp$Y

  # 	Apply the constraints of the general optimization problem
  # 	and compute the outer components
  iter <- 1
  crit <- numeric(n_iter_max)

  a_old <- a
  crit_old <- sum(C * g(cov2(Y, bias = bias)))


  repeat {
    tmp <- sgcca_update(A, a, Y, bias, na.rm, const, J, n, dg, C)
    a <- tmp$a
    Y <- tmp$Y

    # Print out intermediate fit
    crit[iter] <- sum(C * g(cov2(Y, bias = bias)))

    if (verbose) {
      cat(
        " Iter: ", formatC(iter, width = 3, format = "d"),
        " Fit: ", formatC(crit[iter], digits = 8, width = 10, format = "f"),
        " Dif: ", formatC(crit[iter] - crit_old,
          digits = 8, width = 10,
          format = "f"
        ),
        "\n"
      )
    }
    stopping_criteria <- c(
      drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
      abs(crit[iter] - crit_old)
    )

    if (any(stopping_criteria < tol) | (iter > n_iter_max)) {
      break
    }

    crit_old <- crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  for (q in seq_len(J)) {
    if (ctrl & a[[q]][1] < 0) {
      a[[q]] <- -a[[q]]
      Y[, q] <- pm(A[[q]], a[[q]], na.rm = na.rm)
    }
  }

  crit <- crit[which(crit != 0)]

  if (iter > n_iter_max) {
    stop_rgcca(
      "The SGCCA algorithm did not converge after ", n_iter_max,
      " iterations."
    )
  }
  if (verbose) {
    if (iter <= n_iter_max) {
      message(
        "The SGCCA algorithm converged to a stationary point after ",
        iter - 1, " iterations \n"
      )
    }
    plot(crit, xlab = "iteration", ylab = "criteria")
  }

  l2_sat <- sapply(a, function(x) norm(x, "2"))
  if (max(abs(l2_sat - 1)) > tol) {
    for (i in which(abs(l2_sat - 1) > tol)) {
      if (l2_sat[i] < .Machine$double.eps) {
        warning(
          "Norm2 of the block weight vector #",
          i, " is too small :", l2_sat[i]
        )
      } else {
        nMAX <- length(which(a[[i]] != 0))
        warning(
          "The l2 constraint is not saturated for block #", i,
          ". The sparsity parameter has to be in the range [",
          sqrt(nMAX / pjs[i]),
          ", 1] and is equal to ", sparsity[i], "."
        )
      }
    }
  }

  AVE_inner <- sum(C * cor(Y)^2 / 2) / (sum(C) / 2) # AVE inner model

  result <- list(
    Y = Y, a = a, crit = crit,
    AVE_inner = AVE_inner
  )
  return(result)
}
