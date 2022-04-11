#' Projection onto the intersection of the L1 and L2 balls
#' @param x numeric vector to project onto the intersection of the L1 and
#' L2 balls
#' @param a numeric scalar defining the radius of the L1 ball
#' @return \item{k}{number of iteration.}
#' @return \item{lambda}{threshold value for the soft-thresholding operator.}
#' @return \item{l2_sat}{Boolean value indicating whether the
#' l2-constraint are saturated or not.}
#' @return \item{sol}{Numeric vector giving the projection of x onto the
#' intersection of the L1 and L2 balls if the l2-constraint is not saturated.}
#' @importFrom stats median
#' @noRd
proj_l1_l2 <- function(x, a = 1) {

  # Check if constraints are already satisfied
  norm2_x <- norm(x, type = "2")
  if (norm2_x < .Machine$double.eps) {
    return(list(sol = x, l2_sat = FALSE))
  }
  if (sum(abs(x / norm2_x)) <= a) {
    return(list(k = NaN, lambda = 0, l2_sat = TRUE))
  }

  # The desired a_k cannot be null as the constraints are not already satisfied
  # (cf. previous check). So zero values are removed.
  uneq <- x != 0
  p <- abs(x[uneq])

  # Check for multiple maximum
  MAX <- max(p)
  bMAX <- (abs(p - MAX) <= 1e-13)
  nMAX <- sum(bMAX)

  # If there are multiple maximum, the sparse parameter
  # "a" must be >= sqrt(number of max)
  if (a < sqrt(nMAX)) {
    sol <- x * 0
    idx_MAX_val <- which(abs(x) == MAX)
    sol[idx_MAX_val] <- sign(x[idx_MAX_val]) * a / nMAX
    return(list(sol = sol, l2_sat = FALSE))
  }

  # If there are multiple maximum and a = sqrt(number of max),
  # solution is straightforward
  if (a == sqrt(nMAX)) {
    if (nMAX == length(p)) {

      # Case where there is as many maximum as the number of elements of the
      # vector to project "x". Indeed in the case,
      # ||x||_1/||x||_2 = sqrt(number of max)
      lambda <- 0
    } else {

      # With this choice of lambda, the soft-thresholding will set all
      # parameters to 0 except for the elements equals to the maximum value.
      a_2 <- max(p[-which(bMAX)])
      lambda <- (MAX + a_2) / 2
    }
    return(list(k = NaN, lambda = lambda, l2_sat = TRUE))
  }

  # If the vector to project "x" is composed of 2 elements only, as the
  # sparse parameter a <= 1 the desired a_k is a_2 (the lowest element) as
  # psi(a_2) = 1 (cf. theory) and by construction a_3 = 0 and psi(0) >= a
  # (checked previously). Moreover, these 2 elements are necessarily different
  # (cf. conditions above)
  if (length(p) == 2) {
    a_k <- min(p)
    psi_a_k <- 1
    k <- 2
    lambda <- a_k - (a * sqrt((k - psi_a_k^2) / (k - a^2)) - psi_a_k) *
      (sum(p) - k * a_k) / (psi_a_k * (k))
    return(list(k = NaN, lambda = lambda, l2_sat = TRUE))
  }

  # Initialize parameters
  s_1 <- s_2 <- nb <- 0
  repeat {
    N <- length(p)
    if (N == 0) {
      warning("length(p) = 0")
      break
    }

    # Choose next a_k
    if (N %% 2 == 0) {

      # remove max instead of min -> psi not defined for max(abs(x))
      p_reduced <- p[-which.max(p)]

      # Replace by ccaPP::fastMedian or any homemade cpp median function
      a_k <- median(p_reduced)
    } else {

      # Replace by ccaPP::fastMedian or any homemade cpp median function
      a_k <- median(p)
    }

    # Make a partition of list p
    p_inf_ak <- p < a_k
    p_sup_ak <- p > a_k
    p_high <- p[p_inf_ak]
    p_low <- p[p_sup_ak]

    # Evaluation decreasing rank of a_k
    nb_a_k <- sum(p == a_k)
    k <- nb + sum(p_sup_ak) + nb_a_k

    # Compute value of the constraint
    aksq <- a_k^2
    s_low_1 <- sum(p_low) + nb_a_k * a_k

    # NOTE : could create  : ssq   <- function(u) sum(u**2) -> not necessary
    # could use norm(u, type = "2")^2 -> not working when u = as.numeric(0)
    # (mean p_low is empty)
    # When p_low is empty, sum(p_low**2) = 0, which is what is wanted.
    s_low_2 <- sum(p_low**2) + nb_a_k * aksq
    psi_a_k <- (s_1 + s_low_1 - k * a_k) /
      sqrt(s_2 + s_low_2 - 2 * a_k * (s_1 + s_low_1) + k * aksq)

    # Choose partition depending on the constraint
    if (psi_a_k > a) {
      if (length(p_low) == 0) break
      p <- p_low
    } else {
      if (length(p_high) == 0) {
        break
      } else {
        a_k_1 <- max(p_high)
        psi_a_k_1 <- (s_1 + s_low_1 - k * a_k_1) /
          sqrt(s_2 + s_low_2 - 2 * a_k_1 * (s_1 + s_low_1) + k * a_k_1^2)
        if (psi_a_k_1 > a) {
          break
        }
        p <- p_high
        nb <- k
        s_1 <- s_1 + s_low_1
        s_2 <- s_2 + s_low_2
      }
    }
  }

  # Compute lambda
  lambda <- a_k - (a * sqrt((k - psi_a_k^2) / (k - a^2)) - psi_a_k) *
    (s_1 + s_low_1 - k * a_k) / (psi_a_k * (k))
  return(list(k = k, lambda = lambda, l2_sat = TRUE))
}
