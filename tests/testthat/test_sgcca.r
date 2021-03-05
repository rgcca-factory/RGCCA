# '# Test rgcca
#
# '''
#setwd("~/Bureau/RGCCA/tests/testthat")
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(X_agric, X_ind, X_polit);
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);

resSgcca = sgcca(A, C, ncomp=rep(2,3),sparsity= c(1, 1, 1), scheme = "factorial", scale = TRUE,verbose=FALSE)
resRgcca = rgcca(A, connection=C, ncomp=rep(2,3),method="rgcca", scheme = "factorial", scale = TRUE,verbose=FALSE)



resSgcca=sgcca(A,C,sparsity=rep(0.8,3))
head(resSgcca$Y[[3]])
A1=list(pol=A[[3]],agr=A[[1]],ind=A[[2]])
C1=matrix(c(0,1,1,1,0,0,1,0,0),3,3)
resSgcca1=sgcca(A1,C1,sparsity=rep(0.8,3))
head(resSgcca1$Y[[1]])
#
# names(A_final)
# sparsity= c(0.22, 0.08, 0.07, 0.11)
# Cfinal=matrix(0,4,4);Cfinal[1,2:4]=1;Cfinal[2:4,1]=1;
# sparsity2=c(sparsity[2:4],sparsity[1])
# Cfinal2=matrix(0,4,4);Cfinal2[4,1:3]=1;Cfinal2[1:3,4]=1;
# A_final2=list(lipid=A_final[[2]],metabo=A_final[[3]],transcripto=A_final[[4]],cli=A_final[[1]])
# res1=sgcca(A=A_final,C=Cfinal,sparsity=sparsity)
# res2=sgcca(A=A_final2,C=Cfinal2,sparsity=sparsity2)
# head(res1$Y[[1]])
# head(res2$Y[[4]])
#
#
# res1k=sgccak(A=A_final,C=Cfinal,sparsity=sparsity)
# res2k=sgccak(A=A_final2,C=Cfinal2,sparsity=sparsity2)
# head(res1k$Y)
# head(res2k$Y)
# sum(abs(res2k$Y[,4]-res1k$Y[,1]))
#
# # OK les resultats sont les mêmes pour sgccak
# A_final_s=scaling(A_final,scale=TRUE,scale_block=TRUE,bias=FALSE)
# A_final2_s=scaling(A_final2,scale=TRUE,scale_block=TRUE,bias=FALSE)
#
#
# all.equal(A_final_s[[1]],A_final2_s[[4]])
# all.equal(A_final_s[[2]],A_final2_s[[1]])
# all.equal(A_final_s[[3]],A_final2_s[[2]])
# all.equal(A_final_s[[4]],A_final2_s[[3]])
#
# A_final_s=scaling(A_final,scale=TRUE,scale_block=TRUE,bias=FALSE)
# A_final2_s=scaling(A_final2,scale=TRUE,scale_block=TRUE,bias=FALSE)
#
# A_final_s_na=scaling(intersection_list(A_final),scale=TRUE,scale_block=TRUE,bias=FALSE)
# A_final2_s_na=scaling(intersection_list(A_final2),scale=TRUE,scale_block=TRUE,bias=FALSE)
#
# res1k=sgccak(A=A_final_s_na,C=Cfinal,sparsity=sparsity)
# res2k=sgccak(A=A_final2_s_na,C=Cfinal2,sparsity=sparsity2)
#
# head(res1k$Y)
# head(res2k$Y) #pas les bons
# # Without missing values
#
# res1k2=sgccak_2(A=A_final_s_na,C=Cfinal,c1=sparsity)
# res2k2=sgccak_2(A=A_final2_s_na,C=Cfinal2,c1=sparsity2)
#
# head(res1k2$Y)
# head(res2k2$Y)
#
# head(resRgcca$Y[[3]])
# resRgcca1=rgcca(blocks=A1,connection=C1)
# head(resRgcca1$Y[[1]])
