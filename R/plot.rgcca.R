#' Plot for RGCCA
#' 
#' Plot different outputs of the results obtained by a rgcca function
#' @inheritParams plot_var_2D
#' @inheritParams plot_var_1D
#' @inheritParams plot2D
#' @param x A RGCCA object (see \code{\link[RGCCA]{rgcca}})
#' @param ... additional graphical parameters
#' @param type A character among 'ind', 'var', 'both', 'ave', 'cor', 'weight', 
#' 'network', 'factor', 'weight_matrix' (see details).
#' @param text_ind boolean value indicating if rownames are plotted 
#' (default = TRUE). 
#' @param text_var Boolean value indicating if variable names are plotted 
#' (default = TRUE) 
#' @param overlap Boolean value enables avoiding overlapping between labels 
#' (default = FALSE)
#' @param block A vector indicating the blocks to consider.
#' @inheritParams plot_ind
#' @inheritParams plot2D
#' @inheritParams plot_var_2D
#' @inheritParams plot_histogram
#' @details 
#' \itemize{
#' \item "ind" for sample plot. The blocks (block argument) and components 
#' (comp) that will be used on the horizontal and the vertical axis to plot the 
#' individuals: (Y[[block[1]]][, comp[1]], Y[[block[2]]][,comp[2]]). Points can 
#' be colored according to the resp argument. The colors of the points can be 
#' modified with the colors argument.
#' \item  "var" for correlation circle.  
#' first axis, in ordinate, the correlation with the second axis. 
#' \item "both": displays both sample plot and correlation circle (implemented
#' only for one block and at least when two components are asked (ncomp >= 2)
#' \item "ave": displays the average variance explained for each block
#' \item "net": displays the network of connection between blocks (defined by 
#' the connection argument) used in the rgcca() function.
#' \item "cor": barplot of the correlation between variables of 
#' one block (specified by block) and one of its block components (comp). 
#' Variables are sorted in decreasing correlations and only the highest 
#' correlations are displayed. The number of displayed correlations can be set 
#' with n_marks (defaut value = 30).
#' \item "weight": barplot of the block weight vector for one 
#' specific block/component. The weights are sorted from the highest to 
#' the lowest and only the highest are displayed. The number of displayed 
#' weights can be set with n_marks
#' \item "factor": barplot of the block weight factor for one 
#' specific block/component/mode/rank. The weights are sorted from the highest 
#' to the lowest and only the highest are displayed. The number of displayed 
#' weights can be set with n_marks
#' \item "weight_matrix": heatmap of the weight matrix for 3D tensor blocks.
#' }
#' @examples
#' data(Russett)
#' status = colnames(Russett)[9:11][apply(Russett[, 9:11], 1, which.max)]
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
#' X_ind = as.matrix(Russett[,c("gnpr","labo")]);
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
#' A = list(X_agric = X_agric, X_ind = X_ind, X_polit = X_polit);
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
#' fit.rgcca=rgcca(blocks = A, connection = C, 
#'                 tau = rep(1, 3), ncomp = rep(2, 3))
#'                 
#' ###############               
#' # sample plot #
#' ###############
#' 
#' # Defaut call: First component of the last block vs second component of the 
#' # last block
#' plot(fit.rgcca, type="ind", resp = status)
#' 
#' # horizontal axis: First component of the first block
#' # vertical axis: First component of the second block
#' plot(fit.rgcca, type="ind", block = 1:2, comp = 1, resp = status)
#' 
#' # horizontal axis: First component of the first block
#' # vertical axis: Second component of the first block
#' plot(fit.rgcca, type="ind", block = 1, comp = 1:2, resp = status)
#'
#' ######################
#' # Correlation circle #
#' ######################
#' # with superblock
#' fit.mcia = rgcca(blocks=A, scheme = "factorial", ncomp = rep(2, 4), 
#'                  tau = c(1, 1, 1, 0), superblock = TRUE)
#' plot(fit.mcia, type="both", resp = status, overlap = FALSE)
#'                                                      
#' plot(fit.rgcca, type="cor")
#' plot(fit.rgcca, type="weight")
#' plot(fit.rgcca, type="ind")
#' plot(fit.rgcca, type="var")
#' plot(fit.rgcca, type="both")
#' plot(fit.rgcca, type="ave")
#' plot(fit.rgcca, type="network")
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot
#' @export
plot.rgcca=function(x, type = "weight", block = length(x$blocks), comp = 1:2,
                    factor = 1, rank = 1,
                    resp = rep(1, NROW(x$Y[[1]])), remove_var = FALSE, 
                    text_var=TRUE,text_ind=TRUE, response_name = "Response",
                    overlap = TRUE, title = NULL, n_mark = 30, 
                    collapse = FALSE, cex = 1, cex_sub = 12, cex_main = 14, 
                    cex_lab = 12, cex_axis = 10, colors = NULL, ...)
{
    stopifnot(is(x, "rgcca"))
    match.arg(type, c("ind", "var", "both", "ave", "cor", "weight", "network", "factor", "weight_matrix"))
     if(length(comp) == 1){comp = rep(comp,2)}
    compx = comp[1]
    compy = comp[2]

    if(length(block) == 1)
    {
         if(x$call$ncomp[block]<2)
          {
                if(type%in%c("ind", "var", "both"))
                {
                    message("type='ind','var' or 'both' is not available for 
                            ncomp < 2. type was replaced by 'weight'")
                    type="weight"                    
                }

          }
        block=rep(block,2)
  
    }
    i_block=block[1]
    i_block_y=block[2]
    
    for (i in seq(2))
        check_blockx("block", block[1], x$blocks)
      
    if(i_block != i_block_y & is.null(type)){type="weight"}
    if(i_block == i_block_y & is.null(type)){type="both"}
     
    if(type=="both")
    {
        if(is.null(i_block)){i_block=length(x$blocks)}
        p1 <- plot_ind(x, i_block = i_block, i_block_y =i_block_y,
                       compx = compx, compy = compy, cex_sub = cex_sub, 
                       cex_main = cex_main, cex_lab = cex_lab, resp = resp, 
                       response_name = response_name, text = text_ind,
                       title = "Sample space", colors = colors, 
                       no_overlap =!overlap)
        
        p2 <- plot_var_2D(x, i_block = i_block, compx = compx, compy = compy,
                          cex_sub = cex_sub, cex_main = cex_main,
                          cex_lab = cex_lab, remove_var = remove_var,
                          text = text_var, no_overlap = !overlap,
                          title = "Variable correlations", n_mark = n_mark,
                          collapse = collapse, colors = colors)
        
        if(is.null(title)){title = toupper(names(x$blocks)[i_block])}
        p5 <- grid.arrange(p1, p2, nrow=1, ncol=2, top = title)
        return(p5)
    }
    else if(type == "var")
    {
        if(x$call$superblock)
        {
            if(block[1] == length(x$blocks))
            {
                block = length(x$blocks)-1
            }
        }
        if(is.null(title)){title = paste0("Variable correlations: ", 
                                          names(x$blocks)[i_block])}
        
       p5 <- plot_var_2D(x, i_block = i_block, compx = compx, compy = compy,
                         cex_sub = cex_sub, cex_main = cex_main, 
                         cex_lab = cex_lab, remove_var = remove_var,
                         text = text_var, no_overlap =!overlap,
                         title = title, n_mark = n_mark, collapse = collapse,
                         colors = colors)
        return(p5)
    }
    
    else if(type == "ind")
    {
      if(is.null(title))
        {
            if(i_block == i_block_y)
            {
                title= paste0("Sample space: ",names(x$blocks)[i_block])   
            }
            else
            {
                title = "Sample space"
            }
            
        }
        p5<-plot_ind(x, i_block = i_block, i_block_y = i_block_y,
                     compx = compx, compy = compy, cex_sub = cex_sub,
                     cex_main = cex_main, cex_lab = cex_lab, resp = resp,
                     response_name = response_name, text = text_ind,
                     title = title, colors = colors, no_overlap=!overlap)
        return(p5)
     }
    else if(type == "ave")
    {
        if(is.null(title)){title = "Average Variance Explained"}
        p5 <- plot_ave (x, cex = cex, title = title, colors = colors, 
                        cex_main = cex_main, cex_sub = cex_sub)
        return(p5)
    }
    else if(type == "network")
    {
        if(is.null(title)){title=paste0("Common rows between blocks : ",
                                        NROW(x$blocks[[1]]))}
        plot_network ( x, title = title, cex_main = cex_main)
        p5<-NULL
    }
    else if(type == "cor")
    {
        if(is.null(title)){title = paste0("Variable correlations: ", 
                                          names(x$blocks)[i_block])}
        p5=plot_var_1D(x, comp = compx, n_mark = n_mark, type = "cor", 
                       collapse = collapse, title = title, colors = colors, 
                       i_block = i_block, cex_main = cex_main, 
                       cex_sub = cex_sub)
        return(p5)
    }
    else if(type == "weight")
    {
        if(is.null(title)){title= paste0("Variable weights: ",
                                         names(x$blocks)[i_block])}
        
        p5=plot_var_1D(x, comp = compx, n_mark = n_mark, i_block = i_block, 
                       type = "weight", collapse = collapse, title = title, 
                       colors = colors, cex_main = cex_main, cex_sub=cex_sub)
        return(p5)
    }
    else if(type == "factor")
    {
      if (x$call$type != "mgcca") 
        stop_rgcca("Type \"factor\" is only available for MGCCA.")
      nb_dim <- length(dim(x$block[[i_block]]))
      if (factor >= nb_dim) 
        stop_rgcca(paste0("\"factor\" must have a value strictly below ", nb_dim))
      if (rank > x$call$ranks[i_block])
        stop_rgcca(paste0("\"rank\" must not exceed ", x$call$ranks[i_block]))
      if (compx > x$call$ncomp[i_block])
        stop_rgcca(paste0("\"comp\" must not exceed ", x$call$ncomp[i_block]))
      if (nb_dim < 3) type = "weight"
      if(is.null(title)){title= paste0("Variable weights: ",
                                       names(x$blocks)[i_block])}
      
      p5=plot_var_1D(x, comp = compx, n_mark = n_mark, i_block = i_block, 
                     factor = factor, rank = rank,
                     type = type, collapse = collapse, title = title, 
                     colors = colors, cex_main = cex_main, cex_sub=cex_sub)
      return(p5)
    }
    else if(type == "weight_matrix") 
    {
      DIM <- dim(x$block[[i_block]])
      if (length(DIM) != 3) 
        stop_rgcca("Type \"weight_matrix\" is only available for 3D blocks")
      if(is.null(title)){
        title= paste0("Variable matrix weights: ", names(x$blocks)[i_block]) 
      }
      df = expand.grid(
        dimnames(x$blocks[[i_block]])[[2]],
        dimnames(x$blocks[[i_block]])[[3]]
      )
      df$weights = x$a[[i_block]][, compx]
      ggplot(df, aes(df[, 1], df[, 2])) + 
        geom_tile(aes(fill = weights), color = "white") +
        labs(title = title, x = "", y= "", subtitle = print_comp(x, compx, i_block)) +
        #Creating color range
        scale_fill_gradientn(colors=c("skyblue", "yellow", "tomato"), guide="colorbar") +
        #Rotating labels
        theme(axis.text.x = element_text(angle = 270, hjust = 0,vjust=-0.05))
    }

    #invisible(p5)
}
