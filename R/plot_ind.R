#' Plot the two components of a RGCCA
#'
#' Plot the two components of a RGCCA
#'
#' @inheritParams plot2D
#' @param rgcca A list giving the results of a R/SGCCA
#' @param resp A vector of characters corresponding either to a qualitative
#' variable with levels or a continuous variable
#' @param compx An integer giving the index of the analysis component used
#' for the x-axis
#' @param compy An integer giving the index of the analysis component used
#' for the y-axis
#' @param i_block An integer giving the index of a list of blocks
#' @param text A bolean to represent the points with their row names (TRUE)
#' or with circles (FALSE)
#' @param i_block_y An integer giving the index of a list of blocks (another
#' one, different from the one used in i_block)
#' @param reponse_name A character giving the legend title
#' @param no_overlap A boolean to avoid overlap in plotted text
#' @param predicted A list containing as  2nd element a matrix of predicted components
#' @examples
#' coord = lapply(seq(3),
#'    function(x) matrix(runif(15 * 2, min = -1), 15, 2))
#' AVE_X = lapply(seq(3), function(x) runif(2))
#' for (i in 1:length(coord))
#' row.names(coord[[i]]) = seq(15)
#' rgcca_out = list(Y = coord, AVE = list(AVE_X = AVE_X))
#' # Using a superblock
#' resp = as.matrix(rep(LETTERS[seq(3)], each = 5))
#' row.names(resp) = seq(15)
#' rgcca_out$call$type="rgcca"
#' plot_ind(rgcca_out, resp)
#' # Using the first block
#' resp = as.matrix(runif(15, min=-15, max = 15))
#' row.names(resp) = seq(15)
#' rgcca_out$call$type="rgcca"
#' plot_ind(rgcca_out, resp, 1, 2, 1)
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' A = list(X_agric, X_ind, X_polit)
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' result.rgcca = rgcca(A, C, tau = c(1, 1, 1), scheme = "factorial",
#' scale = TRUE,ncomp=rep(2,3))
#' plot_ind(result.rgcca,i_block=1)
#' @export
plot_ind <- function(
    rgcca,
    resp = rep(1, NROW(rgcca$Y[[1]])),
    compx = 1,
    compy = 2,
    i_block = length(rgcca$Y),
    text = TRUE,
    i_block_y = i_block,
    reponse_name = "Response",
    no_overlap = FALSE,
    predicted = NULL,
    title = "Sample space",
    ...){

    if (is.null(i_block_y))
        i_block_y <- i_block

    check_ncol(rgcca$Y, i_block)

    resp <- check_response(resp, rgcca$Y)

    df <- get_comp(
        rgcca = rgcca,
        resp = resp,
        compx = compx,
        compy = compy,
        i_block = i_block,
        i_block_y = i_block_y,
        predicted = predicted
    )
    class(df) <- c(class(df), "d_ind")

    if (!is.null(predicted))
            p <- ggplot(df, aes(df[, 1], df[, 2], color = df$resp))

    else if (length(unique(as.matrix(df$resp))) > 5 && 
            !is.character2(as.vector(df$resp)) ) {

        p <- ggplot(df, aes(df[, 1], df[, 2], color = df$resp))

    }else
        p <- NULL


    p <- plot2D(
            rgcca,
            df,
            title,
            df$resp,
            reponse_name,
            compx,
            compy,
            i_block,
            p,
            text,
            i_block_y,
            no_overlap = no_overlap,
            ...
        )

    # remove legend if missing
    if (length(unique(df$resp)) == 1)
        p <- p + theme(legend.position = "none")

    return(p)
}
