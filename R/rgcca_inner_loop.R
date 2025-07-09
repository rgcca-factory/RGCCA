rgcca_inner_loop <- function(A, C, g, dg, tau = rep(1, length(A)),
                             sparsity = rep(1, length(A)),
                             verbose = FALSE, init = "svd", bias = TRUE,
                             tol = 1e-08, na.rm = TRUE, n_iter_max = 1000,
                             confounders = NULL, penalty_coef = rep(0, length(A)), algo = 1, primal = TRUE) {
  if (!is.numeric(tau)) {
    # From Schafer and Strimmer, 2005
    tau <- vapply(A, tau.estimate, na.rm = na.rm, FUN.VALUE = 1.0)
  }

  # TODO: change this behaviour
  if (any(sparsity == 0)) {
    tau[which(sparsity == 0)] <- 0
    sparsity[which(sparsity == 0)] <- 1
  }

  ### Initialization
  block_objects <- lapply(seq_along(A), function(j) {
    create_block(A[[j]], j, bias, na.rm, tau[j], sparsity[j], tol, confounders[[j]], penalty_coef[j], algo, primal)
  })
  
  block_objects <- lapply(block_objects, block_init, init = init)
  
  Y <- do.call(cbind, lapply(block_objects, "[[", "Y"))
  N <- block_objects[[1]]$N

  iter <- 1
  crit <- NULL
  #crit_RGCCA <- NULL
  #crit_penalty <- NULL
  crit_old <- sum(C * g(crossprod(Y) / N)) 
  crit_blocks <- NULL
  crit_tilde <- NULL
  
  if (!is.null(confounders) && any(penalty_coef != 0)) {
    crit_second_part <- sum(unlist(lapply(seq_along(A), function(j) {
      if (!is.null(confounders[[j]]) && penalty_coef[j] != 0) {
        return(drop(1/N * penalty_coef[j] * t(as.matrix(Y[, j])) %*% confounders[[j]] %*% as.matrix(Y[, j])))
      } else {
        return(0)
      }
    })))
    crit_old <- crit_old - crit_second_part
    crit_blocks <- crit_old
    crit_tilde <- crit_old
  }
  a_old <- lapply(block_objects, "[[", "a")
  if (algo == 5) {z_old <- lapply(block_objects, "[[", "z")}

  repeat {
    for (j in seq_along(A)) {
      # Cat h1
      if (!is.null(confounders) && any(penalty_coef != 0)) {
        h1 <- sum(C * g(crossprod(Y) / N)) - crit_second_part
        #cat("h1 = ", formatC(h1, digits = 8, width = 10, format = "f"), "\n")

        # Cat htilde1
        htilde1 <- sum(C * g(crossprod(Y) / N))  +
          t(2/N * t(block_objects[[j]]$x) %*% Y %*% (C[j, ] * dg(crossprod(Y, Y[, j]) / N))) %*%
          (block_objects[[j]]$a - a_old[[j]]) #- crit_second_part
        cat("htilde1 = ", formatC(htilde1, digits = 8, width = 10, format = "f"), "\n")#, ifelse(h1 == htilde1, "", " NO: h1 != htilde1 "), "\n")
      }
      
      # Compute grad
      grad <- Y %*% (C[j, ] * dg(crossprod(Y, Y[, j]) / N))
      
      # #TEST fct minorante
      # for (i in 1:10) {
      #   w_test <- rnorm(block_objects[[j]]$p)
      #   if (norm(w_test, type = "2") > 1) {w_test <- w_test / norm(w_test, type = "2")}
      #   Y_test <- Y
      #   Y_test[, j] <- block_objects[[j]]$x %*% w_test
      #   f_test <- sum(C * g(crossprod(Y_test) / N))
      #   f_lin <- sum(C * g(crossprod(Y) / N))  + 2/N * t(grad) %*% block_objects[[j]]$x %*% (w_test - block_objects[[j]]$a)
      #   cat("f = ", f_test, " ; f_lin = ", f_lin, ifelse(f_lin - f_test > 1e-10, "NO", ""), "\n")
      # }
      
      block_objects[[j]] <- block_update(block_objects[[j]], grad)
      
      # Cat htilde2
      if (!is.null(confounders) && any(penalty_coef != 0)) {
        crit_second_part <- sum(unlist(lapply(seq_along(A), function(j) {
          if (!is.null(confounders[[j]]) && penalty_coef[j] != 0) {
            return(drop(1/N * penalty_coef[j] * t(block_objects[[j]]$x %*% block_objects[[j]]$a) %*% confounders[[j]] %*% block_objects[[j]]$x %*% block_objects[[j]]$a))
          } else {
            return(0)
          }
        })))}

      if (!is.null(confounders) && any(penalty_coef != 0)) {
        htilde2 <- sum(C * g(crossprod(Y) / N))  +
          t(2/N * t(block_objects[[j]]$x) %*% Y %*% (C[j, ] * dg(crossprod(Y, Y[, j]) / N))) %*%
          (block_objects[[j]]$a - a_old[[j]]) #- crit_second_part
        cat("htilde2 = ", formatC(htilde2, digits = 8, width = 10, format = "f"), ifelse(htilde1 - htilde2 > 1e-10, " NOOOO: htilde1 > htilde2 ", ""), "\n")
      }
      
      Y[, j] <- block_objects[[j]]$Y
      
      # Compute criterion after the block weight vector is updated
      if (!is.null(confounders) && any(penalty_coef != 0)) {
        crit_second_part <- sum(unlist(lapply(seq_along(A), function(k) {
          if (!is.null(confounders[[k]]) && penalty_coef[k] != 0) {
            return(drop(1/N * penalty_coef[k] * t(as.matrix(Y[, k])) %*% confounders[[k]] %*% as.matrix(Y[, k])))
          } else {
            return(0)
          }
        })))
        crit_blocks <- c(crit_blocks, sum(C * g(crossprod(Y) / N))  - crit_second_part)
      } else {
        crit_blocks <- c(crit_blocks, sum(C * g(crossprod(Y) / N)) )
      }

      # Cat h2
      if (!is.null(confounders) && any(penalty_coef != 0)) {
        h2 <- sum(C * g(crossprod(Y) / N))  - crit_second_part
        #cat("h2 = ", formatC(h2, digits = 8, width = 10, format = "f"),  ifelse(htilde2 - h2 > 1e-10, " NOOOOOOOOOOOOOO: htilde2 > h2 ", ""), "\n")
      }
    }
    
    # Print out intermediate fit
    if (!is.null(confounders) && any(penalty_coef != 0)) {
      crit_second_part <- sum(unlist(lapply(seq_along(A), function(j) {
        if (!is.null(confounders[[j]]) && penalty_coef[j] != 0) {
          return(drop(1/N * penalty_coef[j] * t(as.matrix(Y[, j])) %*% confounders[[j]] %*% as.matrix(Y[, j])))
        } else {
          return(0)
        }
      })))
      crit <- c(crit, sum(C * g(crossprod(Y) / N)) - crit_second_part)
      #crit_RGCCA <- c(crit_RGCCA, sum(C * g(crossprod(Y) / N)))
      #crit_penalty <- c(crit_penalty, crit_second_part)
    } else {
      crit <- c(crit, sum(C * g(crossprod(Y) / N)))
    }

    if (verbose) {
      cat(
        " Iter: ", formatC(iter, width = 3, format = "d"),
        " Fit: ", formatC(crit[iter], digits = 8, width = 10, format = "f"),
        " Dif: ", formatC(crit[iter] - crit_old,
                          digits = 8, width = 10, format = "f"),
        #"RGCCA crit: ", formatC(crit_RGCCA[iter], digits = 8, width = 10, format = "f"),
        #"Penalty: ", formatC(crit_penalty[iter], digits = 8, width = 10, format = "f"), 
        "\n"
      )
      
      # for (j in seq_along(A)) {
      #   
      #   cat(
      #     paste("Fit_block_", j, sep = ""), formatC(crit_blocks[length(A) * (iter - 1) + j + 1], digits = 8, width = 10, format = "f"),
      #     paste("Dif_block_", j, sep = ""), formatC(crit_blocks[length(A) * (iter - 1) + j + 1] - crit_blocks[length(A) * (iter - 1) + j], digits = 8, width = 10, format = "f"),
      #     "\n"
      #   )
      # }
    }

    a <- lapply(block_objects, "[[", "a")
    if (algo == 5) {z <- lapply(block_objects, "[[", "z")}
    
    if (algo == 5) {
      stopping_criteria <- c(
        drop(crossprod(unlist(z, FALSE, FALSE) - unlist(z_old, FALSE, FALSE))),
        abs(crit[iter] - crit_old)
      )
    } else {
      stopping_criteria <- c(
        drop(crossprod(unlist(a, FALSE, FALSE) - unlist(a_old, FALSE, FALSE))),
        abs(crit[iter] - crit_old)
        )
    }
    

    if (any(stopping_criteria < tol) || (iter > n_iter_max)) {
      break
    }

    crit_old <- crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  if (iter > n_iter_max) {
    warning(
      "The RGCCA algorithm did not converge after ", n_iter_max,
      " iterations."
    )
  }
  if (verbose) {
    if (iter <= n_iter_max) {
      message(
        "The RGCCA algorithm converged to a stationary point after ",
        iter - 1, " iterations \n"
      )
    }
    plot(crit, xlab = "iteration", ylab = "criteria")
    # plot(crit_blocks, xlab = "iteration", ylab = "criteria", xaxt = "n", col = "green4", pch = 16)
    # axis(side = 1, at = 0:(iter*length(A)), labels = c(0, rep(1:iter, each = length(A))))
    #plot(crit_RGCCA, xlab = "iteration", ylab = "RGCCA criteria")
    #plot(crit_penalty, xlab = "iteration", ylab = "penalty")
  }

  # Post-process the resulting block-weight and block-component vectors
  ctrl <- all(g(-5:5) == g(5:-5))
  block_objects <- lapply(block_objects, block_postprocess, ctrl)
  a <- lapply(block_objects, "[[", "a")
  Y <- do.call(cbind, lapply(block_objects, "[[", "Y"))

  return(list(Y = Y, a = a, crit = crit, tau = tau))
}
