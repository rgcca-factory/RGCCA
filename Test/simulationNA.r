
rm(list=ls())
#------------
# Loading libraries
#-------------------
library(nipals)
library(FactoMineR)
library(RGCCA)
library(missMDA)
library(MASS)
# loading package functions
namesFiles=dir("/home/caroline.peltier/Bureau/RGCCA/R")
# loading functions in R directory
sapply(namesFiles,function(x){source(paste0("/home/caroline.peltier/Bureau/RGCCA/R/",x))})

#--------------------------------
# Creation des jeux de données
#------------------------------------
X1=matrix(rnorm(500),100,5)
X2=matrix(rnorm(1000),100,10)
X1[c(1,4),]=NA
X2[2,]=NA
A=list(X1,X2)
blocknames=c("block1","block2")
nbTestFiles=5
names(A)=blocknames
# Obtenir la probabilité de données manquantes par bloc 
pNA=0.1

# Récupérer le jeu de données complet
A=intersection(A)


# sauver le jdd de reference
wd="/home/caroline.peltier/Bureau"
setwd(wd)
dir.create("Reference")
writeList(A=A,wd=paste0(wd,"/Reference"))
setwd(wd)
# sauver les 20 jeux de données simulés
dir.create("missingValuesSimulation")
setwd(paste0(wd,"/missingValuesSimulation"))
for(i in 1:nbTestFiles)
{
  res1=createNA(A,option="block",pNA=pNA)
  dir.create(as.character(i))
  writeList(res1$dat,wd=as.character(i))
  setwd(paste0(wd,"/missingValuesSimulation"))
}
# read reference data
setwd(paste0(wd,"/Reference"))
refData=readDataset(blocknames)
# Effectuer les analyses de RGCCA et les comparaisons
setwd(paste0(wd,"/missingValuesSimulation"))
an1=analysis(refData,noIntersect=TRUE,C=NULL,tau=NULL,scale=TRUE,nAxe=2,scheme="centroid",sameBlockWeight=TRUE,wd=getwd(),blocknames=blocknames,nbTestFiles=nbTestFiles)

plotAnalysis(an1,ylim=c(0.5,1),output="rv")
plotAnalysis(an1,ylim=c(0.5,1),output="rv")