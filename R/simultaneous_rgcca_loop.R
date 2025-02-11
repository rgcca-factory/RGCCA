#' @importFrom Deriv Deriv
simultaneous_rgcca_loop <- function(blocks, connection = 1 - diag(length(blocks)),
                                    tau = rep(1, length(blocks)),
                                    sparsity = rep(1, length(blocks)),
                                    ncomp = 1,
                                    scheme = "centroid",
                                    init = "svd", bias = TRUE, tol = 1e-08,
                                    verbose = TRUE,
                                    na.rm = TRUE, superblock = FALSE,
                                    response = NULL, disjunction = NULL,
                                    n_iter_max = 1000, comp_orth = TRUE) {
  if (verbose) {
    scheme_str <- ifelse(is(scheme, "function"), "user-defined", scheme)
    cat(
      "Simultaneous computation of the RGCCA block components based on the",
      scheme_str, "scheme \n"
    )
    if (!is.numeric(tau)) {
      cat("Optimal shrinkage intensity parameters are estimated \n")

      # From Schafer and Strimmer, 2005
      tau <- vapply(blocks, tau.estimate, na.rm = na.rm, FUN.VALUE = 1.0)

    }
  }

  if (is.function(scheme)) {
    g <- scheme
  } else {
    switch(scheme,
           "horst" = {
             g <- function(x) x
           },
           "factorial" = {
             g <- function(x) x^2
           },
           "centroid" = {
             g <- function(x) abs(x)
           }
    )
  }

  dg <- Deriv::Deriv(g, env = parent.frame())

  ##### Initialization #####
  J <- length(blocks)
  pjs <- vapply(blocks, NCOL, FUN.VALUE = 1L)
  nb_ind <- NROW(blocks[[1]])

  # Whether primal or dual
  primal_dual <- rep("primal", J)
  primal_dual[which((sparsity == 1) & (nb_ind < pjs))]

  block_objects <- lapply(seq_along(blocks), function(j) {
    is_response <- !is.null(response) && (j == response)
    create_sim_block(blocks[[j]], j, bias, na.rm, tau[j], ncomp[1], is_response)
  })

  block_objects <- lapply(block_objects, block_init, init = init)
  Y <- lapply(block_objects, "[[", "Y")
  N <- block_objects[[1]]$N

  iter <- 1
  crit <- NULL
  crit_old <- sum(vapply(seq(J), function(i) {
    sum(vapply(seq(J), function(j) {
      connection[i, j] * sum(diag(g(crossprod(Y[[i]], Y[[j]]) / N)))
    }, FUN.VALUE = double(1L)))
  }, FUN.VALUE = double(1L)))
  a_old <- lapply(block_objects, "[[", "a")

  ##### Update weights until convergence #####
  repeat {
    for (j in seq_along(blocks)) {
      # Compute grad
      grad <- Reduce("+", lapply(seq(J), function(k) {
        connection[j, k] * Y[[k]] %*% diag(
          dg(diag(crossprod(Y[[j]], Y[[k]]))), nrow = ncol(Y[[k]])
        )
      }))
      block_objects[[j]] <- block_update(block_objects[[j]], grad)
      Y[[j]] <- block_objects[[j]]$Y
    }

    # Print out intermediate fit
    crit <- c(
      crit,
      sum(vapply(seq(J), function(i) {
        sum(vapply(seq(J), function(j) {
          connection[i, j] * sum(diag(g(crossprod(Y[[i]], Y[[j]]) / N)))
        }, FUN.VALUE = double(1L)))
      }, FUN.VALUE = double(1L)))
    )

    if (verbose) {
      cat(
        " Iter: ", formatC(iter, width = 3, format = "d"),
        " Fit: ", formatC(crit[iter], digits = 8, width = 10, format = "f"),
        " Dif: ", formatC(crit[iter] - crit_old,
                          digits = 8, width = 10, format = "f"
        ), "\n"
      )
    }

    a <- lapply(block_objects, "[[", "a")
    stopping_criteria <- c(
      drop(crossprod(unlist(a, FALSE, FALSE) - unlist(a_old, FALSE, FALSE))),
      abs(crit[iter] - crit_old)
    )

    if (any(stopping_criteria < tol) || (iter > n_iter_max)) {
      break
    }

    crit_old <- crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  if (iter > n_iter_max) {
    warning(
      "The simultaneous RGCCA algorithm did not converge after ", n_iter_max,
      " iterations."
    )
  }
  if (verbose) {
    if (iter <= n_iter_max) {
      message(
        "The simultaneous RGCCA algorithm converged to a stationary point after ",
        iter - 1, " iterations \n"
      )
    }
    plot(crit, xlab = "iteration", ylab = "criteria")
  }

  # Post-process the resulting block-weight and block-component vectors
  ctrl <- all(g(-5:5) == g(5:-5))
  block_objects <- lapply(block_objects, block_postprocess, ctrl)
  a <- lapply(block_objects, "[[", "a")
  Y <- lapply(block_objects, "[[", "Y")

  out <- list(
    Y = Y,
    a = a,
    astar = a,
    tau = tau,
    crit = crit, primal_dual = primal_dual
  )

  class(out) <- "grgcca"
  return(out)
}
