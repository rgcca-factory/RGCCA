#' Plot of variables space
#'
#' Correlation circle highlighting the contribution of each variables to the
#' construction of the RGCCA components
#' @inheritParams plot_ind
#' @inheritParams plot2D
#' @param remove_var A bolean to keep only the 100 variables of each
#' component with the biggest correlation#'
#' @param n_mark An integer giving the number of top variables to select
#' @param collapse A boolean to combine the variables of each block as result
#' @param ... Further graphical parameters (see plot2D)
#' @examples
#' setMatrix = function(nrow, ncol, iter = 3) 
#'  lapply(
#'      seq(iter),
#'      function(x) {
#'          y <- matrix(runif(nrow * ncol), nrow, ncol)
#'          rownames(y) <- paste0("S",1:nrow)
#'          return(y)
#'      })
#' blocks = setMatrix(10, 5)
#' blocks[[4]] = Reduce(cbind, blocks)
#' for (i in seq(4)) {
#'  colnames(blocks[[i]]) = paste0( LETTERS[i],
#'  as.character(seq(NCOL(blocks[[i]]))))
#' }
#' coord = setMatrix(10, 2, 4)
#' a = setMatrix(5, 2)
#' a[[4]] = setMatrix(15, 2, 1)[[1]]
#' AVE_X = lapply(seq(4), function(x) runif(2))
#' rgcca_out = list(Y = coord, a = a, AVE = list(AVE_X = AVE_X), call = list(blocks = blocks))
#' names(rgcca_out$a) <- LETTERS[seq(4)] -> names(rgcca_out$call$blocks)
#' rgcca_out$call$type="rgcca"
#' # Using a superblock
#' rgcca_out$call$superblock = TRUE
#' rgcca_out$call$ncomp = rep(3, 4)
#' class(rgcca_out) = "rgcca"
#' plot_var_2D(rgcca_out, 1, 2)
#' # Using the first block
#' plot_var_2D(rgcca_out, 1, 2, 1)
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'  politic = Russett[, 6:11] )
#' rgcca_out = rgcca(blocks,ncomp=2)
#' # Without superblock but with the of all variables to the first block
#' plot_var_2D(rgcca_out, collapse = TRUE)
#' @export
plot_var_2D <- function(
    rgcca_res,
    compx = 1,
    compy = 2,
    i_block = length(rgcca_res$a),
    text = TRUE,
    remove_var = TRUE,
    n_mark = 100,
    collapse = FALSE,
    no_overlap = FALSE,
    title = "Variable correlations",
    resp=NULL,
    colors=NULL,
    ...) {

    x <- y <- NULL
    #if(rgcca_res$call$superblock==FALSE){collapse=FALSE;}
    df <- get_ctr2(
        rgcca_res = rgcca_res,
        compx = compx,
        compy = compy,
        i_block = i_block,
        type = "cor",
        n_mark = n_mark,
        collapse = collapse,
        remove_var = remove_var,
        resp=resp
    )
    class(df) <- c(class(df), "d_var2D")

    if (collapse && rgcca_res$call$superblock) {
        if (i_block == length(rgcca_res$a))
            i_block <- length(rgcca_res$a) - 1
        rgcca_res$a <- rgcca_res$a[-length(rgcca_res$a)]
    }

    if (i_block < length(rgcca_res$a) || tolower(rgcca_res$call$type) == "pca")
        rgcca_res$call$superblock <- FALSE

    # PCA case: remove the superblock in legend
    if (identical(rgcca_res$call$blocks[[1]], rgcca_res$call$blocks[[2]]))
        rgcca_res$call$superblock <- FALSE

    check_ncol(rgcca_res$a, i_block)

    p <- plot2D(
        rgcca_res,
        df,
        title,
        df$resp,
        "Blocks",
        compx,
        compy,
        i_block,
        text = text,
        collapse =  collapse,
        no_overlap = no_overlap,
        colors=colors,
        ...
    )
    
    p <- p+
    geom_path(
        aes(x, y),
        data = plot_circle(),
        col = "grey",
        size =  0.5
    ) +
    geom_path(
        aes(x, y),
        data = plot_circle() / 2,
        col = "grey",
        size = 0.5,
        lty = 2
    )
    
    # remove legend if not on superblock
  
    if ((!rgcca_res$call$superblock || i_block != length(rgcca_res$a)) && !collapse)
    {
        p <- p + theme(legend.position = "none")
    }
    return(p)
}
