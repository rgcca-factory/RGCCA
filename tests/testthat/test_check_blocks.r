#'# imputeEM test

#'''
#'
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab")]);
blocks = list(X_agric,X_ind,X_polit);
blocks2=check_blocks(blocks,add_NAlines=TRUE)
test_that("test_check_unidimensional_blocks",{
    expect_true(rownames(blocks2[[3]])[45]=="Venezuela")
})


X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
colnames(X_ind)=c("gini","labo")
blocks=list(agri=X_agric,ind=X_ind)
blocks2=check_blocks(blocks,add_NAlines=TRUE)


test_that("test_check_same_var_block",{
    expect_true(colnames(blocks2[[1]])[1]=="agri_gini")
})
