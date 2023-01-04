#' Format bootstrap results
#'
#' Raw bootstrap results are stored in a list of n_boot elements where each
#' element of this list is a list of J matrices. Each of these matrices has
#' ncomp[j] columns where ncomp[j] is the number of components for block j.
#' The aim of this function is to reorganize this list into a data.frame.
#' @param W raw bootstrap results
#' @return A data.frame containing the results.
#' @importFrom tidyr pivot_longer
#' @noRd
format_bootstrap_list <- function(W, rgcca_res) {
  ### Combine data.frames for each bootstrap sample, weights and loadings, and
  ### each block.
  Reduce(rbind, lapply(seq_along(W), function(b) {
    Reduce(rbind, lapply(seq_along(W[[b]]), function(i) {
      Reduce(rbind, lapply(seq_along(W[[b]][[i]]), function(j) {
        # For each block, construct a data.frame
        block <- W[[b]][[i]][[j]]
        colnames(block) <- seq_len(NCOL(block))
        df <- data.frame(var = rownames(block), block) %>%
          pivot_longer(!.data$var, names_to = "comp", names_prefix = "X")
        df$block <- names(W[[b]][[i]])[j]
        df$type <- ifelse(i == 1, "weights", "loadings")
        df$boot <- b
        return(df)
      }))
    }))
  }))
}
