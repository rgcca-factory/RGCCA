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



analysis(refData,noIntersect=TRUE,C=NULL,tau=NULL,scale=TRUE,nAxe=2,scheme="centroid",sameBlockWeight=TRUE,filenames)
  

comparison=function(rgcca1,rgcca2,nAxe=1,selec=10,selectPatient=NULL)
{
  diffNorm2=function(vec1,vec2)
  {
    if(!is.null(vec1)&!is.null(vec2))
    {
      c1= cor(vec1,vec2)
      c2= cor(vec1,-vec2)
      return(max(c1,c2))	
    }
    else{ return(NA)}
  }
  if(is.null(selectPatient)){selectPatient=rownames(rgcca1[["Y"]][[1]])}
  J=length(rgcca1$A)
  com=rv=pctBm=rep(NA,J)
  for(i in 1:J)
  {
    refBm=biomarker(resRGCCA=rgcca1,block=i,axes=1,selec=selec)
    com[i]=diffNorm2(rgcca1[["astar"]][[i]][,nAxe],rgcca2[["astar"]][[i]][,nAxe])
    rv[i]=coeffRV(rgcca1[["Y"]][[i]][selectPatient,1:2],rgcca2[["Y"]][[i]][selectPatient,1:2])$rv
    testBm=biomarker(resRGCCA=rgcca2,block=i,axes=1,selec=selec)
    pctBm[i]=sum(names(testBm)%in%names(refBm))/length(refBm)	
  }
  
  return(list(a=com,rv=rv,bm=pctBm))
}


plotEvol=function(listFinale,output="rv",fileName=NULL,ylim=c(0.8,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)
{ #output : "rv", "pct" ou "a"
  #barType="sd" or "stdErr"
  if(is.null(namePlot)){namePlot=output}
  graphics.off()
  if(!is.null(fileName)){png(paste(fileName,".png",sep=""),width=width,height=height)}
  nameData= names(listFinale)
  abscisse=as.numeric(substr(nameData,5,7));names(abscisse)=nameData
  pas=1 
  par(las=1)
  J=length(listFinale[[1]][[1]][[1]][[1]]) #nblock
  if(block=="all"){ split.screen(c(2,2));toPlot=1:J}else{toPlot=block:block}
  namesMethod=names(listFinale[[1]][[1]])
  #colMethod=rainbow(5)[1:length(namesMethod)]
  colMethod=c("cornflowerblue","chocolate1","chartreuse3","red","blueviolet")[1:length(namesMethod)]
  nMeth=0:length(namesMethod)
  names(colMethod)=names(nMeth)=namesMethod
  for(j in toPlot)
  {
    if(block=="all"){screen(j)}
    par(mar=c(5, 4, 4, 2) + 0.1)
    par(mgp=c(3,1,0))
    plot(NULL,main=paste(namePlot,": Block",j),xlim=c(min(abscisse),max(abscisse)+length(namesMethod)-1),ylim=ylim,xlab="Percent of missing values",ylab="Correlation",bty="n")
    axis(side = 1,col="grey",line=0)
    axis(side = 2,col="grey",line=0)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
           "#e9ebec",border="#e9ebec")
    # grid(nx = NULL, ny = NULL, col = "white", lty = "dotted",  lwd = par("lwd"), equilogs = TRUE)
    
    for(da in nameData)
    {
      abline(v=abscisse[da]-1,col="dark grey",lty=2)
      for(rg in namesMethod)
      {
        result=sapply(listFinale[[da]],function(x){return(x[[rg]][[output]][j])})
        moyenne=mean(result)
        if(!barType %in% c("sd","stderr")){ecartType=0}
        if(barType=="sd"){ecartType=sd(result)}
        if(barType=="stderr"){ecartType=sd(result)/sqrt(length(result))}
        
        points(abscisse[da]+pas*nMeth[rg],moyenne,pch=16,col=colMethod[rg])
        segments(abscisse[da]+pas*nMeth[rg],moyenne-ecartType,abscisse[da]+pas*nMeth[rg],moyenne+ecartType,col=colMethod[rg])
      }
    }
    abline(v=abscisse[length(abscisse)]+pas*max(nMeth),col="dark grey",lty=2)
  }
  if(block=="all")
  {
    screen(4)
    legend("center",legend=namesMethod,fill=colMethod,box.lwd=0)
  }
  if(is.numeric(block))
  {
    legend("bottomleft",legend=namesMethod,fill=colMethod,box.lwd=0)
  }
  if(!is.null(fileName)){dev.off()}
}
plotEvols=function(listFinaleNoIntersect,listFinaleIntersect,repName)
{
  setwd("/home/caroline.peltier/Bureau/latex/rgccaEM/")
  dir.create(repName)
  setwd(repName)
  plotEvol(listFinale=listFinaleNoIntersect,output="rv",ylim=c(0.5,1),fileName="rv",block="all",barType="stderr",namePlot="RV")
  plotEvol(listFinale=listFinaleIntersect,output="rv",ylim=c(0.8,1),fileName="rvComplete",block="all",barType="stderr",namePlot="RV on complete patients")
  plotEvol(listFinaleIntersect,output="bm",ylim=c(0.6,1),fileName="topten",block="all",barType="stderr",namePlot="Percent of similar top ten")
  plotEvol(listFinaleIntersect,output="a",ylim=c(0.8,1),fileName="correlation",block="all",barType="stderr",namePlot="Correlations")
}

tau=NULL;C=NULL;scale=TRUE;nAxe=2;noIntersect=TRUE;sameBlockWeight=TRUE;scheme="centroid"

analysis=function(refData,noIntersect=TRUE,C=NULL,tau=NULL,scale=TRUE,nAxe=2,scheme="centroid",sameBlockWeight=TRUE,filenames)
{
  nBlock=length(refData)
  if(is.null(C)) {C=matrix(1,nBlock,nBlock);diag(C)=0;}
  if(is.null(tau)){tau=rep(1,nBlock)}
  listFinale=list()
  library(parallel)
  refRgcca=rgcca(refData[1:nBlock],C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,tau=tau,scale=scale,scheme=scheme,sameBlockWeight=sameBlockWeight,returnA=TRUE)
  for(name in filenames)
  {
    print(name)
    repertoire=paste("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/",name,"/",sep="")
    listFinale[[name]]=mclapply(1:20,mc.cores=7,FUN=function(j)
    {
      print("dataset")
      print(j)
      setwd(paste(repertoire,j,sep=""))
      
      testData=readDataset(c("CLI","MRS","VOL"))
      listRgcca=list()
      #	lAxes=lapply(testData,function(m){resPca=PCA(m,graph=FALSE);eigenValues=resPca$eig[,1];nbAxes=critereCoude(eigenValues,graph=TRUE);return(nbAxes)})	
      #	nbAxesRgcca=as.vector(unlist(lAxes))		
      # listRgcca[["MI-kNN2"]]=MIRGCCA(testData,k=2,niter=5,scale=scale,sameBlockWeight=TRUE,tau,output="weightedMean",scheme=scheme)$rgcca0
      listRgcca[["MI-kNNAll"]]=MIRGCCA(testData,k="all",ni=5,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,output="weightedMean",scheme=scheme,returnA=TRUE,tol=1e-8)$rgcca0
      testDataSB=imputeSB(testData,ncomp=rep(nAxe,nBlock),scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,tol=1e-8,ni=10)
      listRgcca[["EM"]]=rgcca(testDataSB$A,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,scheme=scheme,returnA=TRUE,tol=1e-8)
      listRgcca[["Nipals"]]=rgcca(testData,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,na.rm=TRUE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,scheme=scheme,returnA=TRUE,tol=1e-8)
      if(!noIntersect)
      {
        completeData=intersection(testData)
        listRgcca[["Complete"]]=rgcca(completeData,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,returnA=TRUE,tol=1e-8)
      }
    
      graphics.off()
      if(noIntersect){selectPatient=NULL}else{selectPatient=rownames(listRgcca[["Complete"]][["Y"]][[1]])}
      res=lapply(listRgcca, function(x)
      {
        comp=comparison(rgcca1=refRgcca,rgcca2=x,selectPatient=selectPatient);return(comp)
      }
      )
      
      names(res)=names(listRgcca)
      return(res)
    })
    
  }
  return(listFinale)
}
