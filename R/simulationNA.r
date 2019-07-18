#------------
# Loading libraries
#-------------------
library(nipals)
library(FactoMineR)
library(RGCCA)
library(missMDA)
library(MASS)
# loading package functions
setwd("/home/caroline.peltier/Bureau/EtudeNA")
source("readDataset.r")
source("biomarker.r")
namesFiles=dir("/home/caroline.peltier/Bureau/RGCCA/R")
# loading functions in R directory
sapply(namesFiles,function(x){source(paste0("/home/caroline.peltier/Bureau/RGCCA/R/",x))})

#--------------------------------
# Creation des jeux de données
#------------------------------------
# Obtenir la probabilité de données manquantes par bloc 
pNA=

# Récupérer le jeu de données complet
A=intersection(A)


# sauver le jdd de reference
writeList(A=refData,wd="")

# sauver les 20 jeux de données simulés
dir.create("missingValuesSimulation")
for(i in 1:20)
{
  res1=createNA(A,option="block",pNA=pNA)
  dir.create(as.character(i))
  writeList(res1$dat,repertoire=as.character(i))
}

# Effectuer les analyses de RGCCA et les comparaisons
analysis(refData,noIntersect=TRUE,C=NULL,tau=NULL,scale=TRUE,nAxe=2,scheme="centroid",sameBlockWeight=TRUE,c("missingValuesSimulation"))




