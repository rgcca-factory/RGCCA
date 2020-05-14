# If less than 2 columns, do not run
# x, a list of matrix
# i_block, the position of the tested matrix in the list
# return an error or NULL
check_ncol <- function(x, i_block) {
    if (NROW(x[[i_block]]) < 2) {
        stop_rgcca(
            "This output is available only for more than one variable."
        )
    }
}