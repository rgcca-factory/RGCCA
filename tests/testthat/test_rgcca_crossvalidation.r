#'# rgccacrossvalidation test

#'''
#' data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
 politic = Russett[, 6:11] )
rgcca_out = rgcca.analyze(blocks)
rgcca_crossvalidation(rgcca_out, validation = "kfold", k = 5, n_cores = 1)
 rgcca_crossvalidation(rgcca_out,  validation = "test", n_cores = 1)$scores
rgcca_crossvalidation(rgcca_out, n_cores = 1)

RussettWithNA=Russett
RussettWithNA[1:2,1:3]=NA
RussettWithNA[3,4:5]=NA
RussettWithNA[3,1]=NA
blocksNA = list(agriculture = RussettWithNA[, seq(3)], industry = RussettWithNA[, 4:5],
                politic = RussettWithNA[, 6:11] )

# cross validation
rgcca_out = rgcca.analyze(blocksNA)
#rgcca_crossvalidation(rgcca_out, n_cores = 1)
