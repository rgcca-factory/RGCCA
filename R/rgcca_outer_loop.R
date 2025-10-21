#' @importFrom Deriv Deriv
rgcca_outer_loop <- function(blocks, connection = 1 - diag(length(blocks)),
                             tau = rep(1, length(blocks)),
                             sparsity = rep(1, length(blocks)),
                             ncomp = rep(1, length(blocks)),
                             scheme = "centroid",
                             init = "svd", bias = TRUE, tol = 1e-08,
                             verbose = TRUE,
                             na.rm = TRUE, superblock = FALSE,
                             response = NULL, disjunction = NULL,
                             n_iter_max = 1000, comp_orth = TRUE,
                             rank = 1, mode_orth = 1, separable = TRUE) {
  if (verbose) {
    scheme_str <- ifelse(is(scheme, "function"), "user-defined", scheme)
    cat(
      "Computation of the RGCCA block components based on the",
      scheme_str, "scheme \n"
    )
    if (!is.numeric(tau)) {
      cat("Optimal shrinkage intensity parameters are estimated \n")
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
  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- vapply(blocks, NCOL, FUN.VALUE = 1L)
  nb_ind <- NROW(blocks[[1]])

  crit <- list()
  R <- blocks

  a <- Y <- lambda <- lapply(seq(J), function(b) c())
  factors <- lapply(seq(J), function(b) {
    lapply(seq_along(dim(R[[b]])[-1]), function(m) NULL)
  })

  if (superblock && comp_orth) {
    P <- c()
  } else {
    P <- lapply(seq(J), function(b) c())
  }

  # Save computed shrinkage parameter in a new variable
  computed_tau <- tau
  if (is.vector(tau)) {
    computed_tau <- matrix(
      rep(tau, N + 1),
      nrow = N + 1, J, byrow = TRUE
    )
  }

  if (is.vector(sparsity)) {
    sparsity <- matrix(
      rep(sparsity, N + 1),
      nrow = N + 1, J, byrow = TRUE
    )
  }

  if (is.vector(rank)) {
    rank <- matrix(
      rep(rank, N + 1),
      nrow = N + 1, J, byrow = TRUE
    )
  }

  # Whether primal or dual
  primal_dual <- matrix("primal", nrow = N + 1, ncol = J)
  primal_dual[which((sparsity == 1) & (nb_ind < matrix(
    pjs, nrow = N + 1, ncol = J, byrow = TRUE
  )))] <- "dual"

  ##### Computation of RGCCA components #####
  for (n in seq(N + 1)) {
    if (verbose) {
      cat(paste0(
        "Computation of the RGCCA block components #", n,
        " is under progress...\n"
      ))
    }
    gcca_result <- rgcca_inner_loop(R, connection, g, dg,
                                    tau = computed_tau[n, ],
                                    sparsity = sparsity[n, ],
                                    init = init, bias = bias, tol = tol,
                                    verbose = verbose, na.rm = na.rm,
                                    n_iter_max = n_iter_max,
                                    rank = rank[n, ], mode_orth = mode_orth,
                                    separable = separable
    )

    # Store tau, crit
    computed_tau[n, ] <- gcca_result$tau
    crit[[n]] <- gcca_result$crit

    # Store Y, a, factors and lambda
    a <- lapply(seq(J), function(b) cbind(a[[b]], gcca_result$a[[b]]))
    Y <- lapply(seq(J), function(b) cbind(Y[[b]], gcca_result$Y[, b]))
    factors <- lapply(seq(J), function(b) {
      lapply(seq_along(factors[[b]]), function(m) {
        cbind(factors[[b]][[m]], gcca_result$factors[[b]][[m]])
      })
    })
    lambda <- lapply(
      seq(J), function(b) cbind(lambda[[b]], gcca_result$lambda[[b]])
    )

    # Deflation procedure
    if (n == N + 1) break
    defl_result <- deflate(gcca_result$a, gcca_result$Y, R, P, ndefl, n,
                           superblock, comp_orth, response, na.rm)
    R <- defl_result$R
    P <- defl_result$P
  }

  # If there is a superblock and weight vectors are orthogonal, it is possible
  # to have non meaningful lambda associated to blocks that have been set to
  # zero by the deflation
  if (superblock && !comp_orth) {
    a <- lapply(a, function(x) {
      if (ncol(x) > nrow(x)) {
        x[, seq(nrow(x) + 1, ncol(x))] <- 0
      }
      return(x)
    })
  }

  ##### Generation of the output #####
  if (N == 0) {
    crit <- unlist(crit)
    computed_tau <- as.numeric(computed_tau)
  } else {
    computed_tau <- apply(computed_tau, 2, as.numeric)
  }

  astar <- compute_astar(a, P, superblock, comp_orth, N)

  out <- list(
    Y = Y,
    a = a,
    astar = astar,
    factors = factors,
    lambda = lambda,
    tau = computed_tau,
    crit = crit, primal_dual = primal_dual
  )

  class(out) <- "rgccad"
  return(out)
}
