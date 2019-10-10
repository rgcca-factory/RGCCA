#'Plots multiple imputation for RGCCA
#'
#' This method allows multiple imputation for RGCCA to be represented.
#' @param mi.obj Result of MIRGCCA
#' @param opt.ell "distr" or "ci"
#' @param multiple "ell" or "seg"  representation of multiple points: can be "ell" for ellipses or "seg" for segments
#' @param indnames if TRUE, the names of individuals are displayed
#' @param varnames if TRUE, the names of variables are displayed

#' @param blocks vector of two integers choosing the block(s)represented by the two axes of the block. By default c(1,1) (block 1 for two axes)
#' @param axes  vector of two integers choosing the components in the blocks. By default c(1,2) (first two axes)
#' @param selec number of variabes to be selected in the "fingerprints"
#' @param xlim by default, NULL
#' @param threshod
#' @param cex
#' @param output
#' @param filename
#' @return \item{rgcca0} RGCCA results for the reference dataset
#' @return \item{data} list of imputed data obtained
#' @return \item{rgccaList} list of RGCCA obtained
#' @title plotMIRGCCA: plots the results of MIRGCCA
#' @examples 

plotMIRGCCA=function(mi.obj,opt.ell="distr",multiple="ell",indnames=TRUE,varnames=TRUE,blocks=c(1,1),axes=c(1,2),selec="all",xlim=NULL,threshold=0.5,cex=1,output="R",filename="rgcca.png")
{
  rgcca0=mi.obj$rgcca0
  rgccaList=mi.obj$rgccaList  
  resprocrustes=list()
  nsuj=dim(rgcca0$A[[blocks[1]]])[1]
  niter=length(mi.obj$rgccaList)
  for(i in 1:niter)
  {
    resprocrustes[[i]]=procrustes(X=rgcca0$Y[[blocks[1]]],Y=rgccaList[[i]]$Y[[blocks[1]]],symmetry=TRUE)$Yrot
  }
  corInfo=function(vec,main="",sdV=NULL,cex=0.8)
  {
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",axes=FALSE,main=main)
    text(0.3,seq(0,1,length.out=length(vec)),names(vec),cex=cex)
    if(is.null(sdV))
    {
      text(0.7,seq(0,1,length.out=length(vec)),round(vec,digits=2),col=ifelse(vec>0,"dark green","red"),cex=cex)
    }
    else
    {
      text(0.7,seq(0,1,length.out=length(vec)),paste(round(vec,digits=2),"+/-",sdV,sep=""),cex=cex,col=ifelse(vec>0,"dark green","red"))
      
    }
  }
  graphics.off()
  if(output=="png"){png(filename=filename)}
  split.screen(c(2,2))
  screen(n=1)
  par(mar=c(1,1,3,1))
  if(is.null(xlim))
  {
    minim=min(min(rgcca0$Y[[blocks[1]]][,axes[1]]),min(rgcca0$Y[[blocks[2]]][,axes[2]]))
    maxim=max(max(rgcca0$Y[[blocks[1]]][,axes[1]]),max(rgcca0$Y[[blocks[2]]][,axes[2]]))
  }
  else{minim=xlim[1];maxim=xlim[2]}
  plot(NULL,xlim=c(minim,maxim),ylim=c(minim,maxim),xaxt="n",yaxt="n",main="Sample plot")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#e9ebec",border="#e9ebec")
  grid(nx = NULL, ny = NULL, col = "white", lty = "dotted",  lwd = par("lwd"), equilogs = TRUE)
  nsuj=dim(rgcca0$Y[[1]])[1]
  if(!indnames)
  {
    points(rgcca0$Y[[blocks[1]]][,axes[1]],rgcca0$Y[[blocks[2]]][,axes[2]],pch=16,col=rainbow(nsuj))
  }
  else
  {
    text(rgcca0$Y[[blocks[1]]][,axes[1]],rgcca0$Y[[blocks[2]]][,axes[2]],rownames(rgcca0$Y[[blocks[1]]]),col=rainbow(nsuj),cex=cex)
  }
  abline(v=0,col="grey")
  abline(h=0,col="grey")
  if(multiple=="seg")
  {
    for(i in 1:niter)
    {
      segments(resprocrustes[[i]][,1],resprocrustes[[i]][,2],rgcca0$Y[[blocks[1]]][,1],rgcca0$Y[[blocks[1]]][,2],lty=3,lwd=0.5,col=rep(rainbow(nsuj),niter))
    }
    
  }
  if(multiple=="ell")
  {
    for(k in 1:nsuj)
    {
      tablePoints=NULL
      for(i in 1:niter)
      {
        tablePoints=rbind(tablePoints,c(resprocrustes[[i]][k,1],resprocrustes[[i]][k,2]))
      }
      if(opt.ell=="distr"){radius=1.96}else{radius=1.96/sqrt(niter)}
      car::ellipse(center=c(mean(tablePoints[,1],na.rm=T),mean(tablePoints[,2],na.rm=T)),shape=cov(tablePoints),radius=radius,col=rainbow(nsuj)[k],lwd=0.5,center.pch = FALSE)
    }
  }
  
  plotCircle=function(origine=c(0,0),r=1)
  {
    x=origine[1]+r*cos(seq(0,2*pi,length.out=100))
    y=origine[2]+r*sin(seq(0,2*pi,length.out=100))
    lines(x,y,col="black")
  }
  A0=rgcca0$A
  if(blocks[1]==blocks[2]){superbloc0=A0[[blocks[1]]];varCol="blue"}	
  if(blocks[1]!=blocks[2])
  {	
    superbloc0=cbind(A0[[blocks[1]]],A0[[blocks[2]]]);
    varCol=c(rep("blue",ncol(A0[[blocks[1]]])),rep("red",ncol(A0[[blocks[2]]])))	
  }
  x0=y0=c()
  x=y=matrix(NA,niter,ncol(superbloc0))
  for(i in 1:niter)
  {
    A=mi.obj$rgccaList[[i]]$A
    
    if(blocks[1]==blocks[2]){superbloc=A[[blocks[1]]];varCol="blue"}	
    if(blocks[1]!=blocks[2])
    {	
      superbloc=cbind(A[[blocks[1]]],A[[blocks[2]]]);
      varCol=c(rep("blue",ncol(A[[blocks[1]]])),rep("red",ncol(A[[blocks[2]]])))	
    }
    
    for(k in 1:ncol(superbloc))
    { # correlations des notes obtenues pour la variable avec les scores de Y
      x[i,k]=cor(superbloc[,k],as.matrix(rgcca0$Y[[blocks[1]]][,axes[1]]),use="pairwise.complete.obs")
      y[i,k]=cor(superbloc[,k],as.matrix(rgcca0$Y[[blocks[2]]][,axes[2]]),use="pairwise.complete.obs")
    }
    colnames(x)=colnames(superbloc)
    colnames(y)=colnames(superbloc)
  }
  sdx=apply(x,2,"sd",na.rm=T)
  sdy=apply(y,2,"sd",na.rm=T)
  
  screen(n=2)
  par(mar=c(1,1,3,1))
  indices1=sort(abs(x[1,]),index.return=T,decreasing=TRUE,na.last=TRUE)$ix
  #barplot(rev(x[indices1]),horiz=FALSE,main=paste("Axis 1 (","",")",sep=""))
  
  if(selec=="all"){selectionX=1:dim(x)[2]}
  if(selec=="all"){selectionY=1:dim(y)[2]}
  if(is.numeric(selec)){selectionX=selectionY=1:selec}
  corInfo(rev( (x[1,indices1])[selectionX]),main="Axis 1",sdV=round(rev(sdx[indices1])[selectionX],digits=2),cex=1)
  
  screen(n=3)
  par(mar=c(1,1,3,1))
  indices2=sort(abs(y[1,]),index.return=T,decreasing=TRUE,na.last=TRUE)$ix
  corInfo(rev((y[1,indices2])[selectionY]),main="Axis 2",sdV=round(rev(sdy[indices2])[selectionY],digits=2),cex=1)
  
  screen(n=4)
  par(mar=c(1,1,3,1))
  plot(NULL,xlim=c(-1,1),ylim=c(-1,1),xlab=paste("Comp 1 (",names(rgcca0$A)[1],")",sep=""),ylab=paste("Comp 2(",names(rgcca0$A)[2],")",sep=""),xaxt="n",yaxt="n")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#e9ebec",border="#e9ebec")
  grid(nx = NULL, ny = NULL, col = "white", lty = "dotted",  lwd = par("lwd"), equilogs = TRUE)
  
  plotCircle()
  
  for(k in 1:ncol(superbloc0))
  {
    x0[k]=cor(superbloc0[,k],as.matrix(rgcca0$Y[[blocks[1]]][,axes[1]]),use="pairwise.complete.obs")
    y0[k]=cor(superbloc0[,k],as.matrix(rgcca0$Y[[blocks[2]]][,axes[2]]),use="pairwise.complete.obs")
  }
  names(x0)=colnames(superbloc0)
  names(y0)=colnames(superbloc0)
  if(varnames)
  {
    selectBool=(abs(x0)>threshold)|(abs(y0)>threshold)
    text(x0[selectBool],y0[selectBool],colnames(superbloc0)[selectBool],col=varCol,cex=cex)
  }
  else{points(x0,y0,col=varCol,pch=16)}
  abline(v=0)
  abline(h=0)
  if(output=="png"){dev.off()}
  
}