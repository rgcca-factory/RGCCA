library(parallel)
library(FactoMineR)

# test on random data
#---------------------
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4);
rownames(X1)=rownames(X2)=paste("S",1:70)
A=list(X1,X2);
resultComparison=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2))
resultComparison=whichNAmethod(A,listMethods=c("nipals","mean","em","knn1"),patternNA=rep(0.1,2))

plotAnalysis(resultComparison,output="rmse",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)

resultComparison=whichNAmethod(A,listMethods=c("complete","mean","nipals","knn1","knn4","knnA","knn10","em","superblockEM"),patternNA=rep(0.1,2))
plotAnalysis(resultComparison,output="rv",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)
plotAnalysis(resultComparison,output="rvComplete",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)
plotAnalysis(resultComparison,output="a",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)
plotAnalysis(resultComparison,output="bm",fileName=NULL,ylim=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480)


# Attention ! ! Ne fonctionne pas si trop de valeurs manquantes
resultComparison=whichNAmethod(A,listMethods=c("mean"),patternNA=rep(0.5,2))

#--------
# test on russets
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_agric[c(2,4),]=NA
X_ind[1,]=NA
X_polit[5,1]=NA
A = list(agri=X_agric, ind=X_ind, polit=X_polit)
resultComparison=whichNAmethod(A=A,listMethods=c("complete","mean","nipals","knn1","knn4","knnA","knn10","em","superblockEM"),patternNA=c(0.2,0.2,0.2))

plotAnalysis(resultComparison,output="rv",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480,ylim=NULL)
plotAnalysis(resultComparison,output="rvComplete",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480,ylim=NULL)
plotAnalysis(listFinale=resultComparison,output="a",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480)
plotAnalysis(listFinale=resultComparison,output="bm",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480)



# test on biosca
#---------------
setwd("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/Reference")
refData=readDataset(c("CLI","MRS","VOL"))
resultComparison=whichNAmethod(refData,listMethods=c("complete","mean","nipals","knn1","knn4","knnA","knn10","em","superblockEM"),patternNA=c(0.2,0.2,0.2))

plotAnalysis(resultComparison,output="rv",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480,ylim=NULL)
plotAnalysis(resultComparison,output="rvComplete",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480,ylim=NULL)
plotAnalysis(listFinale=resultComparison,output="a",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480)
plotAnalysis(listFinale=resultComparison,output="bm",fileName=NULL,block="all",barType="sd",namePlot=NULL,width=480,height=480)




