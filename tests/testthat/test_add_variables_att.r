data("Russett")
blocks <- list(agriculture = Russett[, seq(3)],
               industry = Russett[, 4:5],
               politic = Russett[, 6:11])
rgcca_k <- rgcca(lapply(blocks, function(x){return(x[-1,])}), ncomp=1)
rgcca_res <- rgcca(blocks, ncomp=1)

blocks[[1]][, 1] = 0
blocks[[1]][1, 1] = 1

rgcca_k_saved <- rgcca_k
rgcca_k$a <- add_variables_submodel(rgcca_res, rgcca_k$a)
rgcca_k$astar <- add_variables_submodel(rgcca_res, rgcca_k$astar)
res_blocks<- add_variables_data(rgcca_res, rgcca_k$call$blocks)

center_att <- add_variables_attr(rgcca_res,
                                 w=lapply(rgcca_k_saved$call$blocks,
                                          function(i)
                                              attr(i, "scaled:center")),
                                 type = "center")

scale_attr <- add_variables_attr(rgcca_res, lapply(rgcca_k_saved$call$blocks,
                                                   function(i)
                                                       attr(i, "scaled:scale"))
                                 )

for (i in seq(length(rgcca_k$call$blocks)))
{
    attr(rgcca_k$call$blocks[[i]], "scaled:center") <- center_att[[i]]
    attr(rgcca_k$call$blocks[[i]], "scaled:scale") <- scale_attr[[i]]
}
# tests quand il y a une seule variable
data("Russett")
blocks2 <- list(agriculture = Russett[, seq(3)],
                industry = Russett[, 4:5],
                politic = Russett[, 6])

block_i = lapply(blocks2,
                 function(x){
                     if(is.null(dim(x))){
                         return(x[-1])
                     } else {
                         return(x[-1,])
                         }
                     }
                 )

rgcca_k <- rgcca(block_i, ncomp = 1)
rgcca_res <- rgcca(blocks2, ncomp = 1)
res_blocks2<- add_variables_data(rgcca_res, rgcca_k$call$blocks)
