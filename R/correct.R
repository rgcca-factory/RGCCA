# Check that each data block is to all
# See this answer:
# https://math.stackexchange.com/a/551947
# There should work
#' Check that the network is fully connected
#'
#' Given the design matrix, checks that all the blocks are connected between them
#' @param x Design matrix, a symmetric matrix with
#' @return A logical value if it is fully connected or not.
#' @references \url{https://math.stackexchange.com/a/551947}
#' @author Lluís Revilla Sancho
#' @export
#' @examples
#' connections_wrong <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), ncol = 3, nrow = 3)
#' connections_well <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), ncol = 3, nrow = 3)
#' connections_wrong # No connection of block 3 with any other
#' correct(connections_wrong)
#' connections_well
#' correct(connections_well)
#' connections_wrong[3, 3] <- 1 # Connecting block 3 with itself
#' correct(connections_wrong) # Still wrong!
correct <- function(x) {
    if (!isSymmetric(x)) {
        return(FALSE)
    }
    A <- x != 0 # Adjacency
    # Repeat the adjaceny as much as it is needed.
    l <- lapply(seq_len(ncol(A) - 1), function(y){A})
    # Calculate the power (there are more efficient ways but for small matrices it should work)
    red <- Reduce(`%*%`, l, init = A, accumulate = TRUE)
    # Add them up (S)
    final <- Reduce(`+`, red)
    all(final != 0)
}
